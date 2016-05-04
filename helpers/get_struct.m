function mystruct = get_struct(type, nstruct)
% GET_STRUCT retrieve custom data structures.
%   This function is designed as a centralization for the different 
%   complex structures used throughout the program so that they can
%   be edited consistently.
%
%   MYSTRUCT = GET_STRUCT(TYPE, SIZE) returns a matrix of size SIZE 
%   of the custom data structure of type TYPE. SIZE can be multi-
%   dimensional.
%
%   MYSTRUCT = GET_STRUCT(TYPE) returns one structure (SIZE = 1).
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 13.05.2014

  % Set the default size
  if (nargin == 1)
    nstruct = 1;
  end

  % Switch between all the different types
  switch type

    % Structure used to parse the original files (reused in myrecording)
    case 'channel'
      mystruct = struct('color', 1, ...             % Color of the channel (index of 'colors')
                        'compression', 'none', ...  % Compression used for the temporary file
                        'cosmics', false, ...        % Remove the cosmic rays in the image (see imcosmics.m)
                        'detrend', false, ...       % Detrend the image (see imdetrend.m)
                        'file', '', ...             % Path to the original file
                        'fname', '', ...            % Name of the corresponding temporary file
                        'hot_pixels', true, ...     % Remove the hot pixels in the image (see imhotpixels.m)
                        'max', -Inf, ...            % Original maximum value used for rescaling
                        'min', Inf, ...             % Original minimum value used for rescaling
                        'metadata', '', ...         % Recordings metadata
                        'normalize', true, ...      % Normalize the whole stack
                        'type', 'dic');             % Type of channel

    % Structure storing the color codes used in the display
    case 'colors'
      [all_maps, junk, types] = brewermap('list');
      quals = strcmp('Qualitative', types);
      all_maps = [all_maps(~quals); strcat('*', all_maps(~quals))];
      all_maps = {@gray, @hot, @parula, all_maps{:}};
      mystruct = struct('colormaps', {all_maps}, ...     % The colors used for the images
                        'spots', {{'r','k','b','g','k'}}, ...                                   % The colors used for the detections
                        'spots_next', {{'b','r','k','b', 'w'}}, ...                             % The second colors for the detections
                        'status', {{'myg','mky','mkg','mbg','myg'}}, ...                        % The colors for the status of cells
                        'links', {{'y','k','k','y','y'}}, ...                                   % The colors for the links between cells
                        'paths', {all_maps([2:end 1])}, ... % The colors for the paths
                        'text', {{'r','k','b','g','k'}});                                       % The colors for the text in the movies

    % Structure to store detections from segmentations
    case 'detection'
      mystruct = struct('carth', NaN(1, 2), ...     % Cartesian position of the detections (Nx2)
                        'cluster', [], ...          % Temporal cluster containing the detection paths
                        'noise', [], ...            % Parameters of the image noise
                        'properties', []);          % Other properties computed on each detection, depending on the algorithm

    % Structure handling the export options
    case 'exporting'
      mystruct = struct('file_name', '', ...                % The name of the file to create
                        'low_duplicates', true, ...         % Do we use low duplicates paths ?
                        'full_cycles_only', false, ...      % Keep only tracks that both start and end with a division
                        'export_data', true, ...            % Do we export a CSV table of the data ?
                        'data_aligning_type', 'time', ...   % How do we align the paths ? (time/start/end)
                        'export_noise', false, ...          % Export the parameters of the noise in each image
                        'export_movie', false, ...          % Do we export an AVI movie ?
                        'movie_show_index', true, ...       % Do we display the track indexes in the movie ?
                        'movie_show_detection', true, ...   % Do we display the detected radii in the movie ?
                        'movie_show_paths', false, ...      % Do we display the connections of the tracking in the movie ?
                        'movie_show_reconstruction', false);% Do we display the reconstructed image using the detected spots ?

    % The few parameters required to filter the image appropriately
    case 'image_filters'
      mystruct = struct('hot_pixels_threshold', 15, ...     % The threshold used to detect hot pixels using MEAN(pixels) +/- THRESH*STD(pixels)
                        'cosmic_rays_threshold', 15, ...    % The threshold used to detect cosmic rays as being separated from lower values by a gap larger than THRESH*MAD
                        'cosmic_rays_window_size', 10, ...  % The size of the overlapping tiles the algorithm works with
                        'detrend_meshpoints', 32);          % The number of positions over the image used for detrending

    % Structure containing the different parameters required for tracking spots
    case 'image_segmentation'
      mystruct = struct('filter_max_size', 10, ...         % Maximal radius (in um), see filter_spots.m
                        'filter_min_size', 0.75, ...        % Minimal radius (in um), see filter_spots.m
                        'filter_min_intensity', 15, ...    % Minimal intensity (multiplied by the variance of the noise), see filter_spots.m
                        'filter_overlap', 0.75, ...        % Minimal overlap (in percents) for spots to be fused together
                        'detrend_meshpoints', 32, ...      % The number of positions over the image used for detrending
                        'denoise_func', @gaussian_mex, ... % Function used to denoise the image (see imdenoise.m)
                        'denoise_size', -1,          ...   % Parameter used by the denoising function (-1: default value)
                        'denoise_remove_bkg', true, ...    % Removes the background uniform value (estimated using estimate_noise.m) ?
                        'force_estimation', 1, ...         % Will force the estimation of the signal intensity on the raw data
                        'atrous_max_size', 5, ...         % Maximal size of the spots to detect (in um), see imatrous.m
                        'atrous_thresh', 10, ...           % Threshold used to detect a valid spot as brighter than THRESH*MAD
                        'maxima_window', [5 5], ...        % Window size used to detect local maxima
                        'estimate_thresh', 1, ...          % Utilizes only the pixels brighter than this threshold (times noise variance) to estimate the spots
                        'estimate_niter', 15, ...          % Maximal number in the estimation procedure, see estimate_spots.m
                        'estimate_stop', 1e-2, ...         % Stopping criterion for the estimation procedure, see estimate_spots.m
                        'estimate_weight', 0.1, ...        % Convergence weight for the estiamtion procedure, see estimate_spots.m
                        'estimate_fit_position', false);   % Fit also the subpixel position of the spot (slower) ?


    % Structure used to handle the metadata provided by the microscope
    case 'metadata'
      mystruct = struct('acquisition_time', [], ... % Time of frame acquisition
                        'channels', {{}}, ...       % List of acquired channels
                        'channel_index', [], ...    % Channel <-> frame
                        'exposure_time', [], ...    % Frame exposure time
                        'frame_index', [], ...      % Time point <-> frame
                        'plane_index', [], ...      % Plane <-> frame
                        'raw_data', '', ...         % Raw metadata string
                        'z_position', []);          % Z-position of the frame

    % Global structure of a recording and analysis
    case 'myrecording'
      mychannel = get_struct('channel', 0);
      mysegment = get_struct('segmentation', 0);
      mytracks = get_struct('tracking', 0);
      mystruct = struct('channels', mychannel, ...      % Channels of the recording
                        'segmentations', mysegment, ... % Segmentation data
                        'trackings', mytracks, ...      % Tracking data
                        'experiment', '');              % Name of the experiment

    % Global structure of options/parameters for an analysis
    case 'options'
      myfilt = get_struct('image_filters');
      mysegm = get_struct('image_segmentation');
      mytrac = get_struct('spot_tracking');
      mytrkf = get_struct('tracks_filtering');
      mysplt = get_struct('splitting');
      mystruct = struct('analyzed_fields', {{'carth'}}, ... % Fields in mymovie used for the analyis
                        'application', {{''}}, ...          % List of the applications (other than segmentation) that will be performed
                        'config_files', {{}}, ...       % The various configuration files loaded
                        'binning', 1, ...               % Pixel binning used during acquisition
                        'ccd_pixel_size', 6.45, ...       % X-Y size of the pixels in um (of the CCD camera, without magnification)
                        'debug', false, ...                 % Debug mode ON
                        'magnification', 20, ...        % Magnification of the objective of the microscope
                        'normalize', true, ...              % Normalize the results of the analysis onto the reference embryo
                        'spot_tracking', mytrac, ...    % Parameters for tracking the spots
                        'filtering', myfilt, ...        % Parameters for filtering the recordings
                        'tracks_filtering', mytrkf, ... % Parameters for filtering the tracks
                        'pixel_size', -1, ...           % X-Y size of the pixels in um (computed as ccd_pixel_size / magnification)
                        'segmenting', mysegm, ...       % Parameters for segmenting the recordings
                        'segmentation_type', 'dic', ...     % Type of segmentation (dic, markers, all)
                        'split_parameters', mysplt, ... % Parameters for the splitting of touching cells
                        'recompute', true, ...              % Recompute previously computed features (mainly segmentaiton and trackings)
                        'time_interval', 300, ...       % Time interval between frames (in seconds)
                        'verbosity', 2);                % Verbosity level of the analysis (0 null, 1 text only, 2 gui, 3 full with plots)

    % Structure used to segment a channel
    case 'segmentation'
      mydetection = get_struct('detection',0);
      mystruct = struct('denoise', true, ...           % Denoise the segmentation (see imdenoise) ?
                        'detrend', false, ...          % Detrend the segmentation (see imdetrend.m) ?
                        'filter_spots', true, ...      % Filter the spots (see filter_spots.m) ?
                        'detections', mydetection, ... % The structure used to store the resulting detections
                        'type', {{}});                 % The type of segmentation

    % Structure containing the parameters used for splitting touching cells
    case 'splitting'
      mystruct = struct('angle_thresh', 0.23, ...
                        'max_area_diff', 2, ...
                        'max_distance', 14.4, ...
                        'max_overlap', 0.37, ...
                        'max_ratio', 0.52, ...
                        'max_score', 0.047);

    % Structure containing the different parameters required for tracking spots
    case 'spot_tracking'
      mystruct = struct('spot_max_speed', 0.025, ...          % Maximal speed of displacement of a spot (in um/s)
                        'allow_branching_gap', false, ...     % Allows merging/splitting to occur at the same time as gap closing ?
                        'min_section_length', 5, ...          % The minimum number of contiguous frames a path "section" should last to be kept for merging/spliting/gap closing
                        'bridging_max_gap', 3, ...            % Considered number of frames for the gap closing algorithm (see track_spots.m)
                        'max_intensity_ratio', Inf, ...       % Defines an upper bound to the allowed signal ratios (see track_spots.m)
                        'bridging_max_dist', Inf, ...         % Maximal distance throughout the gap (see track_spots.m)
                        'bridging_function', @bridging_cost_sparse_mex, ...   % Function used to measure the gap-closing weight
                        'joining_function', @joining_cost_sparse_mex, ...     % Function used to measure the joinging weight
                        'splitting_function', @splitting_cost_sparse_mex, ... % Function used to measure the splitting weight
                        'linking_function', @linking_cost_sparse_mex); ...    % Function used to measure the frame-to-frame linking

    % The trackings as stored after segmentation
    case 'tracking'
      mydetection = get_struct('detection',0);
      mystruct = struct('reestimate_spots', true, ...      % Do we reestimate the newly interpolated spots ?
                        'filtered', mydetection, ...       % The structure used to store the detections after filtering
                        'detections', mydetection);        % The structure used to store the resulting detections

    % The options for filtering tracks
    case 'tracks_filtering'
      mystruct = struct('interpolate', true, ...        % Interpolate the position of cells to fill gaps ?
                        'max_zip_length', 3, ...        % A "zip" is defined as a cell splitting and merging back together. This defines the maximum number of frames this separation can last to be closed
                        'min_path_length', 5);          % The minimum number of frames a path should last to be kept

    % If the required type of structure has not been implemented, return an empty one
    otherwise
      mystruct = struct();
  end

  % Compute the pixel size
  mystruct = set_pixel_size(mystruct);

  % Repeat the structure to fit the size (nstruct can be multi-dimensional)
  mystruct = repmat(mystruct, nstruct);

  return;
end
