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
    thresh = mean(atrous(:)) + noise_thresh * std(atrous(:));
    % And get the list of candidates
    [cand_x, cand_y] = find(atrous > thresh);
    nspots = length(cand_x);

    % Prepare the parameters for the optimization
    opts = optimset('Display','off', 'Algorithm', 'levenberg-marquardt');
    
    % Initialize the intermediate variables used to store the detections
    results = NaN(nspots, 6, nspots);
    init_conds = NaN(nspots, 7);

    % We then loop over every candidate spot to refine our search
    for j = 1:nspots
      % Get a first estimate of the position and size of the detected spot
      [gauss_params] = estimate_spot(img, [cand_x(j), cand_y(j)], spot_max_size, noise_thresh);

      % If it was not significant, we just skip it
      if (isempty(gauss_params))
        continue;
      end

      % Otherwse we store our initial guess for later comparison
      init_conds(j, 1:end-2) = gauss_params;
      
      % Check whether there was another initial condition close to ours, as similar initial conditions would converge
      % to the same optimization, we can then jsut skip this "costly" step
      indx = find((init_conds(1:j-1, 1) == gauss_params(1)) & (init_conds(1:j-1, 2) == gauss_params(2)), 1);
      dists = hypot(init_conds(1:j-1, 1) - gauss_params(1), init_conds(1:j-1, 2) - gauss_params(2));

      % We use a fixed threshold to define the distance to our initial condition
      indx = find(dists < dist_thresh, 1);

      % If no such "close" initial condition exists, we optimize our new candidate
      if (isempty(indx))

        % We use MATLAB non-linear solver with the gaussian-fitting error function
        [best] = lsqnonlin(@fit_gaussian, gauss_params, [], [], opts);

        % If the final fit has no amplitude (4th parameter), this is not a valid spot
        if (best(4) == 0)
          % Reset the initial conditions not to misslead other candidates, and skip 
          init_conds(j, :) = NaN;
          continue;
        end

        % Otherwise check whether there is a "close" optimized solution
        dists = hypot(results(1:j-1, 1, 1) - best(1), results(1:j-1, 2, 1) - best(2));

        % We merge overlapping spots
        indx = find(dists <= (results(1:j-1, 3, 1) + best(3)), 1);

      % Otherwise we simply copy the optimized result from the table
      else
        sub_indx = init_conds(indx, end);
        indx = init_conds(indx, end-1);
        best = results(indx, 1:end-1, sub_indx);
      end

      % If we have not found a "close" optimized soliton, we store ours as a new one
      if (isempty(indx))
        % Store it in the current row in the first place
        % Adding at the end the "a trous" detection score
        results(j,:,1) = [best atrous(cand_x(j), cand_y(j))];

        % Store the reference indexes used for reusing the solution
        init_conds(j, end-1:end) = [j 1];

      % Otherwise we simply add it in the same row, in the next available place
      else
        % Find the next available free place
        sub_indx = find(isnan(results(indx,1,:)),1);

        % Store the parameters and the indexes
        results(indx,:,sub_indx) = [best atrous(cand_x(j), cand_y(j))];
        init_conds(j, end-1:end) = [indx sub_indx];
      end
    end

    % Average the different optimization results, this is why we copy the optimized results
    % from close initial conditions, so that each dot "weights" accordingly in the final position
    results = mymean(results, 3);

    % Remove NaN solutions
    if (~isempty(results))
      spot = results(~isnan(results(:, 1)), :);

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
  function err = fit_gaussian(params)
  % The function returns the error between the current estimation of the spot and the actual
  % image as a simple sum of absolute difference. 

    % Retrieve the current value of the parameters
    tmp_pos = params(1:2);
    tmp_sigma = params(3);
    tmp_ampl = params(4);
    tmp_bkg = params(5);

    % Get the closest pixel as the center of the mask is on the central pixel (help GaussMask2D)
    pix_pos = round(tmp_pos);

    % Get the size of the window around the candidate spot (less noise and more flat so that it's 
    % easier to fit au gaussian)
    wsize = ceil(1.5*tmp_sigma);

    % Extract a window of size wsize x wsize around pix_pos in the full image
    window = get_window(img, pix_pos, wsize);

    % Compute the gaussian spot as estimated using the provided parameters
    gauss = GaussMask2D(tmp_sigma, size(window), pix_pos - tmp_pos);
    % The spot need to be rescaled properly
    gauss = gauss*tmp_ampl + tmp_bkg;

    % Compute the absolute difference between the two. An addition penality term is added as
    % the 'levenberg-marquardt' optimization algorithm does not handle boundaries on the parameters
    % while negative values do not make any sense in our problem.
    err = sum(sum(abs(gauss - window))) + exp(-10*sum(params(params < 0))) - 1;

    return;
  end
end

function [gauss_params] = estimate_spot(img, estim_pos, wsize, noise_thresh)
% This function performs the initial esitmation of the size and shape of the spot

  % Extract a window around the candidate spot
  spots = get_window(img, estim_pos, wsize);

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
  window = zeros(2*wsize + 1);

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
