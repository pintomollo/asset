function spots = detect_spots(imgs, opts)
% SPOT_DETECT Generic algorithm that detects spots in biological images using the
% IMATROU function as a basis. Using these candidates, it keeps only the strongest
% detections and fits a gaussian signal around it and merges overlapping spots.
%
%   SPOTS = DETECT_SPOTS(IMG, OPTS) returns a list of detect spots in IMG using the
%   parameters provided by OPTS. It uses the parameters in the field 'spot_tracking'
%   (get_struct('spot_trackings')) along with the 'pixel_size' (help set_pixel_size).
%   SPOTS is a Nx6 matrix where each row has the structure [x y s A b p]:
%     - x,y   the subpixel mean of the gaussian
%     - s     the sigma of the gaussian
%     - A     the amplitude of the gaussian
%     - b     the background signal of the gaussian
%     - p     the detection score as returned by IMATROU
%
%   SPOTS = DETECT_SPOTS(IMG) uses the default values as provided by get_struct('ASSET').
%
%   SPOTS = DETECT_SPOTS(IMGS, ...) returns a cell-vector containing the detections 
%   for each plane contained in IMGS.
%
%   SPOTS = DETECT_SPOTS(MYMOVIE, ...) performs the detection in every plane of the
%   recording contained in the structure provided by get_struct('mymovie').
%
% References:
%   [1] Olivo-Marin, J.-C. Extraction of spots in biological images using 
%       multiscale products. Pattern Recognit. 35, 1989-1996 (2002).
%   
%   [2] Jaensch, S. et al. Automated tracking and analysis of centrosomes in early 
%       Caenorhabditis elegans embryos. Bioinformatics. (2010)
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 10.12.2010

  global myspots;

  % Input checking
  if (nargin < 2)
    opts = get_struct('ASSET');
  end

  % Adapting to possibly different types of images
  % If it's a string, we assume it's a filename, if it's a structure we assume it's a 'mymovie'
  if (ischar(imgs) |isstruct(imgs))
    % In this case, we store the filename and load the number of frames from it
    fname = imgs;
    nframes = size_data(fname);

  % Otherwise we initialize the filename and treat imgs as a stack
  else
    fname = '';
    nframes = size(imgs, 3);
  end

  % Initialize the output variable
  spots = cell(nframes, 1);

  % Get the maximum spot size (in um), used only to reduce the amount of computation
  spot_max_size = opts.spot_tracking.max_size;

  % If the size of the puxels has been set, we can compute the actual spot size
  if (opts.pixel_size > 0)
    
    % The maximal size of the spot in pixels
    spot_max_size = ceil(spot_max_size / opts.pixel_size);
  else
    
    % Try to compute the pixel size
    opts = set_pixel_size(opts);

    % If it worked, we can now compute the size in pixels
    if (opts.pixel_size > 0)
      spot_max_size = ceil(spot_max_size / opts.pixel_size);
    else
      % Otherwise we take it as such
      spot_max_size = ceil(spot_max_size);
    end
  end

  % Get the two parameters used in the detection
  % The first one is used to remove noise from the signal
  noise_thresh = opts.spot_tracking.noise_thresh;
  % The second one is to fuse close spots
  if (opts.pixel_size > 0)
    dist_thresh = opts.spot_tracking.fusion_thresh / opts.pixel_size;
  end

  % We iterate over the frames
  for i = 1:nframes
    
    % If we have a file, we load the image
    if (~isempty(fname))
      img = imnorm(double(load_data(fname, i)));

    % Otherwise we get the plane
    else
      img = imgs(:, :, i);
    end
 
    % Initial filtering the image as proposed by [2] to remove noise and increase signal
    % as this ressembles a deconvolution step
    img = imfilter(img, fspecial('gaussian', 7, 0.6), 'symmetric');

    % Performs the actual spot detection "a trous" algorithm [1]
    atrous = imatrou(img, noise_thresh, spot_max_size);

    % Keeps only the detections that are above the noise_threshold
    thresh = median(atrous(:)) + noise_thresh * std(atrous(:));

    bw = imfill(atrous > thresh, 'holes');
    bw = bwmorph(bw, 'shrink', Inf);

    % And get the list of candidates
    [cand_x, cand_y] = find(bw);
    nspots = length(cand_x);

    % Prepare the parameters for the optimization
    opt = optimset('Display','off', 'Algorithm', 'levenberg-marquardt');
    %{
    opt = cmaes('defaults');
    opt.MaxFunEvals = 200;
    opt.TolFun = 1e-5;
    opt.SaveFilename = '';
    opt.SaveVariables = 'off';
    opt.EvalParallel = 'no';
    opt.LogPlot = 0;
    opt.DispModulo = Inf;
    %}

    if (opts.spot_tracking.recursive)
      orig_img = img;
      size_img = size(img);
    end

    % Initialize the intermediate variables used to store the detections
    %results = NaN(nspots, 6, nspots);

    %% Use the last one as an index for mymean
    results = NaN(nspots, 6 + 1);
    init_conds = NaN(nspots, 4);

    % We then loop over every candidate spot to refine our search
    for j = 1:nspots

      % Get a first estimate of the position and size of the detected spot
      [gauss_params] = estimate_spot(img, [cand_x(j), cand_y(j)], spot_max_size, noise_thresh);

      % If it was not significant, we just skip it
      if (isempty(gauss_params))
        continue;
      end

      % Otherwise we store our initial guess for later comparison
      init_conds(j, 1:end-1) = gauss_params;
      
      % Check whether there was another initial condition close to ours, as similar initial conditions would converge
      % to the same optimization, we can then jsut skip this "costly" step
      
      %indx = find((init_conds(1:j-1, 1) == gauss_params(1)) & (init_conds(1:j-1, 2) == gauss_params(2)), 1);
      [dists] = hypot(init_conds(1:j-1, 1) - gauss_params(1), init_conds(1:j-1, 2) - gauss_params(2));

      
      % We use a fixed threshold to define the distance to our initial condition
      indx = find(dists < dist_thresh, 1); 

      % If no such "close" initial condition exists, we optimize our new candidate
      if (isempty(indx) | isnan(init_conds(indx, end)))

        %rescaling = ones(size(gauss_params));
        %rescaling = 10.^(floor(log10(gauss_params)));
        %rescaling(gauss_params == 0) = 1;

        %wsize = ceil(3*gauss_params(3));
        %window = get_window(img, round(gauss_params(1:2)), wsize);

        params_range = 6*gauss_params([3 3 3])/(2*pi);
        init_params = [gauss_params(1:2) gauss_params(3)];

        [junk, myindx] = min(hypot(myspots(:, 1) - gauss_params(1), myspots(:, 2) - gauss_params(2)));

        %params_range = [gauss_params([3 3]) * 2 spot_max_size]/pi;
        %init_params = [gauss_params(1:2) spot_max_size/2];

        first = true;
        % We use MATLAB non-linear solver with the gaussian-fitting error function
        %[bps] = lsqnonlin(@fit_gaussian, [0 0 tan((gauss_params(3) - init_params(3)) / params_range(3))], [], [], opt);
        [bps] = lsqnonlin(@fit_gaussian, [0 0 0], [], [], opt);
        [err, best] = fit_gaussian(bps);
        %keyboard

        %gauss_params = best;
        %wsize = ceil(2*gauss_params(3));
        %window = get_window(img, round(gauss_params(1:2)), wsize);
        %[best] = lsqnonlin(@fit_gaussian, gauss_params, [], [], opt);
        %[best] = cmaes(@fit_gaussian, gauss_params(:), 0.5, opt);
        %best = best(1:2) .* rescaling(1:2);
        %keyboard
        
        % If the final fit has no amplitude (4th parameter), this is not a valid spot
        if (best(4) == 0)
          % Reset the initial conditions not to misslead other candidates, and skip 
          init_conds(j, :) = NaN;
          %display('kept')
          %continue;
        end

        % Otherwise check whether there is a "close" optimized solution
        dists = hypot(results(1:j-1, 1) - best(1), results(1:j-1, 2) - best(2));

        % We merge overlapping spots
        %%%indx = find(dists <= (results(1:j-1, 3) + best(3)), 1);
        indx = [];

      % Otherwise we simply copy the optimized result from the table
      else
        %sub_indx = init_conds(indx, end);
        indx = init_conds(indx, end);
        best = results(indx, 1:end-2);
      end

      % If we have not found a "close" optimized soliton, we store ours as a new one
      if (isempty(indx))
        % Store it in the current row in the first place
        % Adding at the end the "a trous" detection score
        results(j,:) = [best atrous(cand_x(j), cand_y(j)) j];

        % Store the reference indexes used for reusing the solution
        init_conds(j, end) = j;

        %%%%%%%%% CHECK IF REMOVING NOW IS BETTER !!!!!!!!


        if (opts.spot_tracking.recursive & best(4) > 0)
          img = img - GaussMask2D(best(3), size_img, best(1:2), 0, 1)*best(4);

          imagesc(img)
          drawnow
        end

      % Otherwise we simply add it in the same row, in the next available place
      else
        display('found')
        % Find the next available free place
        %sub_indx = find(isnan(results(indx,1,:)),1);

        % Store the parameters and the indexes
        results(j,:) = [best atrous(cand_x(j), cand_y(j)) indx];
        %init_conds(j, :) = NaN;
        init_conds(j, end) = indx;
      end
    end

    % Average the different optimization results, this is why we copy the optimized results
    % from close initial conditions, so that each dot "weights" accordingly in the final position
    
    results = mymean(results(:, 1:end-1), 1, results(:, end));
    %temp_res = mymean(results(:, 1:end-1, :), 3);
    %results = cat(2, temp_res, mysum(results(:, 6, :), 3));
    
    if (opts.spot_tracking.recursive)

      opts.spot_tracking.recursive = false;
      new_spots = detect_spots(img, opts);
      results = [results; new_spots];

      keyboard
    end

    % Remove NaN solutions
    if (~isempty(results))
      spot = results(~any(isnan(results), 2), :);

      % If there is at least one spot, we sort them according to their "a trous" score
      if (~isempty(spot))
        [junk, indx] = sort(spot(:, end), 1, 'descend');
        spot = spot(indx,:);

        % In case there is only one frame, we return the spots directly
        if (nframes == 1)
          spots = spot;

        % Otherwise we create a cell vector
        else
          spots{i} = spot;
        end
      end
    end

    if (iscell(spots) & isempty(spots{i}))
      spots{i} = NaN(0, 6);
    elseif (isempty(spots))
      spots = NaN(0, 6);
    end
  end


  return;

  % The error function used to fit the gaussian spot estimation, I used a nested function to
  % avoid having to pass the whole image as a parameter.
  function [err, params] = fit_gaussian(params)
  % The function returns the error between the current estimation of the spot and the actual
  % image as a simple sum of absolute difference. 

    params = atan(params) .* params_range + init_params;

    % Retrieve the current value of the parameters
    %params(1:2) = atan(params(1:2)) * pos_range + init_pos;
    %params(3) = exp(params(3));

    tmp_pos = params(1:2);
    tmp_sigma = params(3);

    if (tmp_sigma <= 0.5)
      err = Inf;

      return;
    end

    %tmp_sigma = gauss_params(3);
    %tmp_ampl = params(4);
    %tmp_bkg = params(5);
    %tmp_pos = params(1:2) .* rescaling(1:2);
    %tmp_sigma = gauss_params(3) * rescaling(3) / 2;
    %tmp_sigma = params(3) * rescaling(3);
    %tmp_ampl = params(4) * rescaling(4);
    %tmp_bkg = params(5) * rescaling(5);

    % Get the closest pixel as the center of the mask is on the central pixel (help GaussMask2D)
    pix_pos = round(tmp_pos);
    %pix_pos = round(gauss_params(1:2));

    % Get the size of the window around the candidate spot (less noise and more flat so that it's 
    % easier to fit au gaussian)
    wsize = ceil(1.5*tmp_sigma);

    % Extract a window of size wsize x wsize around pix_pos in the full image
    window = get_window(img, pix_pos, wsize);

    % Compute the gaussian spot as estimated using the provided parameters
    gauss = GaussMask2D(tmp_sigma, size(window), tmp_pos - pix_pos);

    frac = 0.98;
    strength = gauss / (2*pi*tmp_sigma^2);
    strength = frac*(strength) + (1-frac)/numel(gauss);

    frac = 0.5;
    goods = (gauss > exp(-frac));

    %scatter(gauss(:), window(:))
    %keyboard

    if (sum(goods) > 10)
      c = [ones(sum(goods(:)), 1), gauss(goods)] \ window(goods);
    else
      c = [ones(numel(gauss), 1), gauss(:)] \ window(:);
    end

    %if (any(c < 0))

    %  err = Inf;
    %  params = [params 0 0];

    %  return;
    %end
    c(c < 0) = 0;

    % The spot need to be rescaled properly
    %gauss = gauss*tmp_ampl + tmp_bkg;
    gauss = gauss*c(2) + c(1);

    % Compute the absolute difference between the two. An addition penality term is added as
    % the 'levenberg-marquardt' optimization algorithm does not handle boundaries on the parameters
    % while negative values do not make any sense in our problem.
    %err = sum(sum(abs(gauss - window))) + exp(-10*sum(params(params < 0))) - 1;
    %smallers = (window < gauss);
    %err = (1.25*sum(gauss(smallers) - window(smallers)) + sum(window(~smallers) - gauss(~smallers)) + exp(-10*sum(params(params < 0))) - 1) / numel(window);
    %err = (sum(sum(abs(gauss - window))) + exp(-10*sum(params(params < 0))) - 1) / numel(window);
    err = sum((abs(gauss(:) - window(:)) .* strength(:))); % + exp(-10*sum(params(params < 0))) - 1;
    err = mymean(abs(gauss(:) - window(:))); % .* strength(:))); % + exp(-10*sum(params(params < 0))) - 1;

    %smallers = (window < gauss);
    %err = (20*sum((gauss(smallers) - window(smallers)) .* strength(smallers)) + sum((window(~smallers) - gauss(~smallers)) .* strength(~smallers))); % + exp(-10*sum(params(params < 0))) - 1) / numel(window);
    params = [params c([2 1]).'];

    if (false)
    %figure;
    a = tmp_pos-pix_pos+wsize+1-tmp_sigma/2;
    subplot(1,2,2)
    hold off;
    %imagesc(strength)
    %imagesc(img)
    imagesc(gauss);
    hold on;
    rectangle('Position', [tmp_pos([2 1])-tmp_sigma/2 tmp_sigma tmp_sigma], 'Curvature', [1 1], 'EdgeColor', 'k');
    subplot(1,2,1);
    hold off;
    %imagesc(window)
    imagesc(gauss - window)
    hold on;
    rectangle('Position', [a([2 1]) tmp_sigma tmp_sigma], 'Curvature', [1 1], 'EdgeColor', 'k');
    title(num2str(err))
    drawnow
    pause(0.1 + 0.5*first)
    first = false;
    end

    if (false)
    a = tmp_pos-pix_pos+wsize+1-tmp_sigma/2;
    b = myspots(myindx, 1:2)-pix_pos+wsize+1-myspots(myindx, 3)/2;
    hold off;
    imagesc(window)
    %imagesc(gauss - window)
    hold on;
    rectangle('Position', [a([2 1]) tmp_sigma tmp_sigma], 'Curvature', [1 1], 'EdgeColor', 'k');
    rectangle('Position', [b([2 1]) myspots(myindx, [3 3])], 'Curvature', [1 1], 'EdgeColor', 'w');
    title(num2str(abs(myspots(myindx, 1:3) - params(1:3)) ./ params_range));
    drawnow
    pause(0.1)
    end

    %keyboard

    return;
  end
end

function [gauss_params] = estimate_spot(img, estim_pos, wsize, noise_thresh)
% This function performs the initial estimation of the size and shape of the spot

  % Extract a window around the candidate spot
  spots = get_window(img, estim_pos, wsize);

  mask = GaussMask2D(wsize/2, size(spots)); 
  spots = spots .* mask;

  indx = [-wsize:wsize];
  [X,Y] = meshgrid(indx.', indx);
  rads = hypot(X(:), Y(:));
  signs = sign(X(:));
  signs(signs == 0) = 1;

  x = rads .* signs;
  y = spots(:);

  y = y - min(y);
  %y = y / sum(y);

  %sigma1 = sqrt(sum(y .* (x.^2)))

  y = y / max(y);
  sigma = mymean(real(sqrt(-x.^2 ./ (2*log(y)))));

  center_x = wsize+1;
  center_y = wsize+1;

  if (false)
  keyboard
  test = elliptic_coordinate(spots, [wsize; wsize]/2, [wsize; wsize], 0);

  pos = -wsize:wsize;
  x = sum(spots, 1);
  y = sum(spots, 2).';

  x = x - min(x);
  y = y - min(y);

  x = x / max(x);
  y = y / max(y);

  [junk, center_x] = max(x);
  [junk, center_y] = max(y);

  sigma_x = (sum(x > exp(-0.5)*max(x)) - 1) / 2;
  sigma_y = (sum(y > exp(-0.5)*max(y)) - 1) / 2;

  sigma = (sigma_x + sigma_y)/2;
  end

  if (false)
    hold off;
    imagesc(spots)
    hold on
    rectangle('Position', [([center_x center_y])-sigma/2 sigma sigma], 'Curvature', [1 1], 'EdgeColor', 'k');
    rectangle('Position', [([center_x center_y])-sigma1/2 sigma1 sigma1], 'Curvature', [1 1], 'EdgeColor', 'r');
    drawnow
    pause(1)
  end

  center_x = center_x + estim_pos(2) - (wsize+1);
  center_y = center_y + estim_pos(1) - (wsize+1);

  gauss_params = [center_y center_x sigma];

  return;

  % Segment the image using the standard deviation and the provided threshold
  bw = (spots > mymean(spots(:)) + noise_thresh*mad(spots(:)));

  % Get the center and other useful information from this segmentation
  props = regionprops(bw, spots, 'Centroid', 'EquivDiameter', 'MaxIntensity', 'MinIntensity');

  % If we did not get anything, just return an empty vector
  if (length(props) == 0)
    gauss_params = [];
    return;

  % If there is only one candidate, we choose it
  elseif (length(props) == 1)
    indx = 1;

  % Otherwise, we'll choose the candidate that is closer to the initial estim_pos
  else
    % Initialize the index and distance
    indx = 1;
    d = hypot(props(1).Centroid(1) - wsize - 1, props(1).Centroid(2) - wsize - 1);

    % Loop over all the candidates
    for i = 2:length(props)

      % Get the distance to the current candidate
      tmp = hypot(props(i).Centroid(1) - wsize - 1, props(i).Centroid(2) - wsize - 1);

      % If it's closer than our actual best, replace it
      if (tmp < d)
        d = tmp;
        indx = i;
      end
    end
  end

  % Get the estimated center of the spot
  pos = props(indx).Centroid;
  pos = pos(1, [2 1]);
  pos = pos - wsize - 1 + estim_pos;

  % Get the estimated standard deviation of the gaussian spot
  sigma = props(indx).EquivDiameter / 3;

  % The intensity of the background
  bkg = props(indx).MinIntensity;

  % The amplitude of the spot
  ampl = props(indx).MaxIntensity - bkg;

  % Organize the parameters correctly
  gauss_params = [pos, sigma, ampl, bkg];
  
  return;
end

function window = get_window(img, pos, wsize)
% This function extracts a window of size wsize in img around pos.

  % Get the size of the full image
  [h,w] = size(img);

  % Initialize the window
  window = NaN(2*wsize + 1);

  % Create an index vector to compute which pixels lie outside the image 
  indx = [1:2*wsize+1];

  % First check along the x coordinate
  indxx = indx - wsize - 1 + pos(1);
  % We keep only the valid indexes
  okx = (indxx > 0 & indxx <= h);
  % Then along the y
  indxy = indx - wsize - 1 + pos(2);
  oky = (indxy > 0 & indxy <= w);

  % Assign the pixels from the image to the window. All the outside ones will be NaN
  window(indx(okx), indx(oky)) = img(indxx(okx), indxy(oky));

  return;
end
