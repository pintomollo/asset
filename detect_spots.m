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
  if (ischar(imgs))
    % In this case, we store the filename and load the number of frames from it
    fname = imgs;
    [nframes, size_img] = size_data(fname);
  
    % Initialize the output variable
    spots = cell(nframes, 1);

  elseif (isstruct(imgs) & isfield(imgs, 'experiment'))
    
    mymovie = imgs;
    fname = mymovie.data.fname;
    [nframes, size_img] = size_data(fname);

    if (isfield(mymovie.data, 'spots') & ~isempty(mymovie.data.spots))
      spots = mymovie.data.spots;
    else
      spots = get_struct('ruffles', [1,nframes]);
    end

  % Otherwise we initialize the filename and treat imgs as a stack
  else
    fname = '';

    size_img = size(imgs);
    nframes = size_img(3);

    size_img = size_img(1:2);

    % Initialize the output variable
    spots = cell(nframes, 1);
  end


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

  cpb = ConsoleProgressBar();
  cpb.setLeftMargin(4);   % progress bar left margin
  cpb.setTopMargin(1);    % rows margin
  cpb.setLength(100);      % progress bar length: [.....]
  cpb.setMinimum(0);
  cpb.setMaximum(nframes);

  cpb.setElapsedTimeVisible(1);
  cpb.setRemainedTimeVisible(1);
  cpb.setElapsedTimePosition('left');
  cpb.setRemainedTimePosition('right');

  cpb.start();

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
      img = double(load_data(fname, i));

    % Otherwise we get the plane
    else
      img = imgs(:, :, i);
    end
 
    params = estimate_noise(img);

    % Initial filtering the image as proposed by [2] to remove noise and increase signal
    % as this ressembles a deconvolution step
    img = imfilter(img, fspecial('gaussian', 7, 0.6), 'symmetric');

    % Performs the actual spot detection "a trous" algorithm [1]
    atrous = imatrou(img, spot_max_size);

    [atAvg, atStd] = localAvgStd2D(atrous, 3*spot_max_size + 1);
    thresh = (atrous > (atAvg + noise_thresh * atStd) & img > (params(1) + noise_thresh*params(2)));
    thresh = bwmorph(thresh, 'clean');
    bw = locmax2d(img .* thresh, ceil(spot_max_size([1 1])/3));

    % And get the list of candidates
    [cand_y, cand_x] = find(bw);
    nspots = length(cand_x);

    results = NaN(nspots, 5);

    % We then loop over every candidate spot to refine our search
    for j = 1:nspots

      % Get a first estimate of the position and size of the detected spot
      [gauss_params] = estimate_spot(img, [cand_x(j), cand_y(j)], spot_max_size, params);

      % If it was not significant, we just skip it
      if (isempty(gauss_params))
        continue;
      end

      % Otherwise we store our initial guess
      results(j, :) = [gauss_params atrous(cand_y(j), cand_x(j))];

      % And remove it from the image
      img = img - GaussMask2D(gauss_params(3), size_img, gauss_params([2 1]), 0, 1)*gauss_params(4);
    end

    % Remove NaN solutions
    if (~isempty(results))
      spot = results(~any(isnan(results), 2), :);

      % If there is at least one spot, we sort them according to their "a trous" score
      if (~isempty(spot))
        [junk, indx] = sort(spot(:, end), 1, 'descend');
        spot = spot(indx,:);

        % Store the results in the corresponding structure if possible
        if (isstruct(spots))

          spots(i).carth = spot(:, 1:2);
          spots(i).properties = spot(:, 3:end);

        % In case there is only one frame, we return the spots directly
        elseif (nframes == 1)
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

    text = sprintf('Progress: %d/%d', i, nframes);
    cpb.setValue(i);
    cpb.setText(text);
  end

  cpb.stop();

  if (isstruct(spots))
    mymovie.data.spots = spots;
    spots = mymovie;
  end

  return;
end

function [gauss_params] = estimate_spot(img, estim_pos, wsize, noise_params)
% This function performs the initial estimation of the size and shape of the spot

  gauss_params = [];

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

  y = y - noise_params(1);

  ampl = max(y);

  if (ampl < 0 | ~isfinite(ampl))
    return;
  end

  y = y / ampl;

  goods = (y > 0 & y < 1);
  if (~any(goods))
    return;
  end

  sigma = mymean(real(sqrt(-x(goods).^2 ./ (2*log(y(goods))))));

  gauss_params = [estim_pos sigma ampl];

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
