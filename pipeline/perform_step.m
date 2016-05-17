function [varargout] = perform_step(cast_step, segment_type, varargin)
% PERFORM_STEP performs the steps specific to the selected segmentation approach,
% thus providing a flexible and modular approach for CAST to analyze various types
% of recordings.
%
%   [...] = PERFORM_STEP(CAST_STEP, SEGMENT_TYPE, ...) performs the CAST_STEP of the
%   pipeline for the SEGMENT_TYPE of approach. The additional inputs as well as the
%   output variables are dependent on which CAST_STEP is to be performed.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 31.03.2015

  % Which step are we performing ?
  switch cast_step

    % Detect embryos in an image
    case 'embryo'

      % We'll need an image and some options
      img = double(varargin{1});
      opts = varargin{2};

      % Now determine which approach to use and which threshold to set
      use_edges = false;
      switch segment_type
        case 'dic'
          thresh = 0.5;
          use_edges = true;
        case 'eggshell'
          img = max(img(:)) - img;
          noise_params = estimate_noise(img);
          thresh = noise_params(1) + noise_params(2);
        case 'fluorescence'
          noise_params = estimate_noise(img);
          noise_thresh = min(round(range(img(:)) / (15*noise_params(2))) + 1, 10);
          thresh = noise_params(1) + noise_thresh*noise_params(2);
        otherwise
          varargout = {[], []};
          return;
      end

      % First estimate where there are "stuff" in the image
      estimates = detect_embryos(img, use_edges, ...
                  opts.split_parameters.egg_min_size / opts.pixel_size, thresh, ...
                  opts.split_parameters.border_thresh);

      % Then split this estimate into elliptical embryos
      embryos = split_cells(estimates, opts.split_parameters.angle_thresh, ...
                opts.split_parameters.max_ratio, ...
                opts.split_parameters.max_distance / opts.pixel_size, ...
                opts.split_parameters.max_score, opts.split_parameters.max_overlap);

      % The ouput variables
      varargout = {embryos, estimates};

    % We will here estimate a set of detections
    case 'estimation'

      % For that we need an image, spots and the options
      img = varargin{1};
      spots = varargin{2};
      opts = varargin{3};

      % Check if we are forcing the estimation
      force_estim = false;
      if (length(varargin) > 3)
        force_estim = varargin{4};
      end

      % Now determine which approach to use
      switch segment_type

        % Gaussian spots
        case 'multiscale_gaussian_spots'
          if (force_estim)
            spots = estimate_spots(img, spots, opts.segmenting.atrous_max_size/(2*opts.pixel_size), []);
          else
            spots = estimate_spots(img, spots, opts.segmenting.atrous_max_size/(2*opts.pixel_size), ...
                               opts.segmenting.estimate_thresh, ...
                               opts.segmenting.estimate_niter, ...
                               opts.segmenting.estimate_stop, ...
                               opts.segmenting.estimate_weight, ...
                               opts.segmenting.estimate_fit_position);
          end

        % Estimation windows
        case 'rectangular_local_maxima'
          spots = estimate_window(img, spots, opts.segmenting.maxima_window);

        % We do not know what to do here...
        otherwise
          spots = [];
      end

      % Return the estimated spots
      varargout = {spots};

    % Here we segment an image
    case 'segmentation'

      % The specific inputs
      img = varargin{1};
      opts = varargin{2};
      noise = varargin{3};

      % The various approaches
      switch segment_type
        case 'multiscale_gaussian_spots'
          window_size = opts.segmenting.atrous_max_size / opts.pixel_size;
          min_sigma = (opts.segmenting.filter_min_size / opts.pixel_size)^2;
          avg_spot = (((sqrt(2*pi*min_sigma))/(2*window_size)) * ...
                        erf(window_size / sqrt(2*min_sigma)))^2;

          min_intens = opts.segmenting.filter_min_intensity*noise(2) + noise(1);

          spots = detect_spots(img, opts.segmenting.atrous_thresh, ...
                               window_size, min_intens);
        case 'rectangular_local_maxima'
          spots = detect_maxima(img, opts.segmenting.maxima_window);
        otherwise
          spots = [];
      end

      % The ouput variable
      varargout = {spots};

    % Here we compute the total intensity of the estimated spots
    case 'intensity'

      % The specific input
      spots = varargin{1};

      % The various approaches
      switch segment_type
        case 'multiscale_gaussian_spots'
          spots_intens = intensity_gaussians(spots);
        case 'rectangular_local_maxima'
          spots_intens = intensity_windows(spots);
        otherwise
          spots_intens = NaN(size(spots, 1), 1);
      end

      % The output
      varargout = {spots_intens};

    % Here we filter the spots according to their properties, so we need to sort
    % out the various extrema values
    case 'filtering'

      % The specific inputs
      spots = varargin{1};
      opts = varargin{2};
      noise = varargin{3};

      % Check if we are forcing the estimation
      re_estim = false;
      if (length(varargin) > 3)
        re_estim = varargin{4};
      end

      % The various approaches are used to determine the proper inputs for the
      % actual filtering function
      switch segment_type
        case 'multiscale_gaussian_spots'
          spots_intens = intensity_gaussians(spots);
          extrema = [opts.segmenting.filter_min_size / opts.pixel_size, ...
                     opts.segmenting.filter_min_intensity*noise(2); ...
                     opts.segmenting.filter_max_size / opts.pixel_size, Inf];
          fusion = @fuse_gaussians;

          % Now lower intensity threshold for re-estimation
          if (re_estim)
            extrema(1,2) = 0;
          end
        case 'rectangular_local_maxima'
          spots_intens = intensity_windows(spots);
          extrema = [opts.segmenting.filter_min_size / opts.pixel_size, ...
                     opts.segmenting.filter_min_intensity*noise(2); ...
                     opts.segmenting.filter_max_size / opts.pixel_size, Inf];
          if (size(extrema, 2) < 3)
            extrema = extrema(:, [1 1:end]);
          end
          fusion = @fuse_windows;

          % Now lower intensity threshold for re-estimation
          if (re_estim)
            extrema(1,3) = 0;
          end
        otherwise
          spots = [];
          spots_intens = [];
          fusion = [];
          extrema = [];
      end

      % Filter the spots using the intensities, extrema and fusion specific to
      % each segmentation approach
      [filtered, goods] = filter_spots(spots, spots_intens, fusion, extrema, ...
                                       opts.segmenting.filter_overlap);

      % The ouput variables
      varargout = {filtered, goods};

    % Here we reconstruct an image based on the detections
    case 'reconstructing'

      % The specific inputs
      orig_img = varargin{1};
      spots = varargin{2};

      % The various approaches
      switch segment_type
        case 'multiscale_gaussian_spots'
          % We take advantage of the GaussMask2D library function for that !
          draw = @(params,ssize)(GaussMask2D(params(3), ssize, params([2 1]), 0, 1) * params(4));
        case 'rectangular_local_maxima'
          draw = @draw_window;
        otherwise
          draw = [];
          spots = [];
      end

      % Reconstruct the image
      img = reconstruct_detection(orig_img, real(spots), draw);

      % The ouput variable
      varargout = {img};

    % Draw the detected spots
    case 'plotting'

      % The specific inputs
      handle = varargin{1};
      spots = varargin{2};
      colors = varargin{3};

      % The various approaches
      switch segment_type
        case 'multiscale_gaussian_spots'
          hgroup = plot_gaussians(handle, spots, colors);
        case 'rectangular_local_maxima'
          hgroup = plot_windows(handle, spots, colors);
        otherwise

          % Otherwise, we need to create an empty group of handles to be consistent
          if (strncmp(get(handle, 'Type'), 'hggroup',7))
            hgroup = handle;
            delete(get(hgroup, 'Children'));
          else
            hgroup = hggroup('Parent', handle);
          end
      end

      % The output
      varargout = {hgroup};

    % Define the values for exporting the data of the spots
    case 'exporting'

      % The inputs
      opts = varargin{1};
      int_scale = varargin{2};

      % The various approaches
      switch segment_type
        case 'multiscale_gaussian_spots'
          % The names of the properties
          colname = {'status', 'x_coord_um', 'y_coord_um', 'sigma_um', 'amplitude_int', 'interpolated'};
          % The factors for the various conversions
          rescale_factor = [1 ([1 1 1] * opts.pixel_size) int_scale 1];
        case 'rectangular_local_maxima'
          colname = {'status', 'x_coord_um', 'y_coord_um', 'width_um','height_um', 'mean_int', 'standard_deviation_int', 'interpolated'};
          rescale_factor = [1 ([1 1 1 1] * opts.pixel_size) ([1 1] *int_scale) 1];
        otherwise
          colname = {};
          rescale_factor = [];
      end

      % The ouputs
      varargout = {colname, rescale_factor};

    otherwise
      varargout = {};
  end

  return;
end
