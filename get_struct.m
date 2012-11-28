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
% 09.12.2010

  % Set the default size
  if (nargin == 1)
    nstruct = 1;
  end

  % Switch between all the different types
  switch type

    % The general structure that contains all the information required for
    % ASSET to analyze the data. Usually this structure is called 'opts' throughout the code.
    case 'ASSET'
      mystruct = struct('analyzed_fields', {{'carth'}}, ... % Fields in mymovie used for the analyis
                        'application', {{''}}, ...          % List of the applications (other than segmentation) that will be performed
                        'auto_save', true, ...              % Automatically save the intermediate results
                        'binning', 1, ...                   % Pixel binning used during acquisition
                        'ccd_pixel_size', 6.45, ...         % X-Y size of the pixels in µm (of the CCD camera, without magnification)
                        'compression', 'none', ...           % Compression used for the data files (prompted is empty)
                        'compute_probabilities', false, ... % Compute the posterior probability
                        'config_file', '', ...              % Name of the configuration file that will be loaded 
                        'crop_export', false, ...           % Crop the images when exporting the results
                        'crop_size', 2.2, ...               % Crop size is the axes_length*crop_size
                        'debug', false, ...                 % Debug mode ON
                        'do_ml', 'none', ...                % Machine learning (ML) is performed
                        'dp_method', 'double', ...          % Dynamic programming method used (see dynamic_programming.m)
                        'export_movie', false, ...          % Export the results of the analysis
                        'file_regexpr', '(.+[-_])?(.*?[-_]?)(\.[\w\.]+)', ... % The regular expression representing the files
                        'filters', get_struct('channel_filter'), ...
                        'follow_periphery', true, ...       % Follow only the periphery of the trackings (no invaginations)
                        'force_circularity', true, ...      % Enforces that the last row of the DP conincides with the first one
                        'magnification', 63, ...            % Magnification of the objective of the microscope
                        'max_export', 1, ...                % Number of exported frames (>1 is the absolute number of frames, <=1 is the fraction)
                        'max_frames', 1, ...                % Number of analyzed frames (value as for max_export)
                        'measure_performances', false, ...  % Measures the error of the segmentation with respect to the manual trackings
                        'merge_input_files', true, ...      % Uses the Bio-formats file merge when importing movie files
                        'ml_type', 'cortex', ...            % Field on which ML is performed (use cell array for more than one)
                        'nbins', 36, ...                    % Number of bins used to measure the performance
                        'normalize', true, ...              % Normalize the results of the analysis onto the reference embryo
                        'overwrite', true, ...              % Overwrite the previous data by saving in the same MAT-file
                        'parse_export', 'normal', ...       % How export is performed (normal or random)
                        'parse_frames', 'normal', ...       % Order of the frames for the segmentation (normal or random)
                        'pixel_size', 0, ...                % X-Y size of the pixels in ï¿½m (computed as ccd_pixel_size / magnification)
                        'quantification', get_struct('quantification'), ... % Parameters of the quantification
                        'recompute', false, ...             % Recompute previously computed features (mainly segmentaiton and trackings)
                        'segment', true, ...                % Perform the segmentation (useful when combine with recompute)
                        'segmentation_parameters', get_struct('segmentations'), ... % Parameters of the segmentation
                        'segmentation_type', 'dic', ...     % Type of segmentation (dic, markers, all)
                        'spot_tracking', get_struct('spot_tracking'), ... % Parameters of the spot tracking algorithm
                        'split_parameters', get_struct('splitting'), ... % Parameters for the splitting of touching cells
                        'temperatures', get_struct('temperatures'), ... % Parameters of the posterior decoding
                        'trackings', '', ...                % List of tracking files
                        'uuid', 0 , ...                     % Universal Unique IDentifier (used in ML to identify processes) 
                        'verbosity', 2, ...                 % Verbosity level (0 null, 1 text only, 2 gui, 3 full with plots)
                        'warp_type', 'radial');             % Warp type used to normalize the embryo (see carth2normalized.m)

      % Compute the pixel size based on the default values. This needs to be re-done
      % in case one value is changed.
      mystruct = set_pixel_size(mystruct);

    % Structure used to parse the original files (reused to create the fields of mymovie)
    case 'channel'
      mystruct = struct('axes_length', zeros(2, 0), ...
                        'centers', zeros(2, 0), ...
                        'color', ones(1,3), ...             % Color of the channel (RGB)
                        'compression', 'none', ...           % Compression used for the temporary file
                        'cortex', [], ...
                        'detrend', false, ...               % Detrend the image (see imdetrend.m)
                        'eggshell', [], ...
                        'file', '', ...                     % Path to the original file
                        'fname', '', ...                    % Name of the corresponding temporary file
                        'hot_pixels', true, ...             % Remove the hot pixels in the image (see imhotpixels.m)
                        'max', -Inf, ...                    % Original maximum value used for rescaling
                        'min', Inf, ...                     % Original minimum value used for rescaling
                        'neighbors', [], ...
                        'planes',1,...                      % Number of planes in a z-stack 
                        'orientations', zeros(1, 0), ...    
                        'type', 'dic', ...                  % Type of channel (dic, eggshell, cortex, data)
                        'timing', get_struct('timing'), ... % Timing of the cell cycle
                        'update', zeros(2, 0));

    case 'channel_filter'
      mystruct = struct('channel', 'data', ...
                        'applied', false, ...
                        'filter', '', ...
                        'params', {{}});

    % Structure used to store the parameters of the correction function (see duplicate_segmentation.m)
    % The correction function is : F(i) = a + b*I(i+s) + c*R(I)
    case 'conversion'
      mystruct = struct('bkg', 0.0424, ...                  % Value of the background shift ('a')
                        'factor', -0.0440, ...              % Value of the intensity factor ('b')
                        'range', -0.0176, ...               % Value of the range factor ('c')
                        'safety', 0, ...                    % Value used to expand the eggshell to avoid excluding part of the membrane signal
                        'shift', 5.1836);                   % Value of the shift between the intensity and the correction ('s')

    % Parameters used to compute the data part of the DP scoring function (see dynamic_programming.m)
    case 'data_parameters'
      mystruct = struct('filt', [], ...                     % Filter applied to the image
                        'path', [], ...                     % Path of interest (usually the eggshell)
                        'alpha', 0, ...                     % Parameters that can be used in the scoring function (7 provided)
                        'beta', 0, ...                      %  |
                        'gamma', 0, ...                     %  |
                        'delta', 0, ...                     %  |
                        'epsilon', 0, ...                   %  |
                        'zeta', 0, ...                      %  |
                        'eta', 0);                          %  |

    % Structure used to store the results of the segmentations (both eggshell and cortex) (see segment_movie.m)
    case {'eggshell', 'cortex'}
      mystruct = struct('carth', [], ...                    % Nx2 matrix of points representing the contour in Cartesian coordinates
                        'estim', [], ...                    % Initial segmentation (without DP, as used in the paper)
                        'temperatures', [], ...             % Temperatures found to compute the posterior probability
                        'thickness', 0, ...                 % Thickness of the eggshell in ER
                        'warped', []);                      % Normalized position

    % Structure used to handle tracking files (see import_trackings.m)
    case 'file'
      mystruct = struct('fname','', ...                     % Path of the original .shapes file (see load_shapes.m)
                        'groups',{{}}, ...                  % Name of the shape groups found in the file (see load_shapes.m)
                        'shapes',[]);                       % Shapes extracted from the file

    case 'fitting'
      mystruct = struct('type', 'simul', ...
                        'init_noise', 0.2, ...
                        'data_noise', 0.2, ...
                        'max_iter', 10000, ...
                        'nfits', 3, ...
                        'fit_full', true, ...
                        'fit_flow', false, ...
                        'parameter_set', 2, ...
                        'scale_type', 'normalize', ...
                        'step_size', 0.01, ...
                        'ground_truth', [], ...
                        'x_pos', [], ...
                        't_pos', [], ...
                        'temperature', 1, ...
                        'fraction', 0.775, ...
                        'aligning_type', 'domain');

    case 'flow'
      mystruct = struct('distance', [], ...
                        'index', [], ...
                        'position', [], ...
                        'speed', []);

    % Structure used to handle the metadata provided by the microscope
    case 'metadata'
      mystruct = struct('acquisition_time', [], ...
                        'channels', {{}}, ...
                        'channel_index', [], ...
                        'exposure_time', [], ...
                        'frame_index', [], ...
                        'plane_index', [], ...
                        'z_position', []);


    % Parameters used to perform machine learning (see find_parameters.m)
    case 'ml_params'                              
      mystruct = struct('config', {{}}, ...                 % Parameters of the optimization
                        'evolution', {{}}, ...              % "Path" of the optimization
                        'goal', [], ...
                        'initial_condition', [], ...        % Starting position for the optimization
                        'ml_type', '', ...                  % Type of segmentation which is optimized (eggshell, cortex)
                        'params', [], ...                   % Current value of the optimized parameters
                        'score', Inf);                      % Score of the current iteration

    case 'modeling'
      mystruct = struct('nparticles', 200, ...
                        'npopulations', 2, ...
                        'dimensionality', 1, ...
                        'x_step', 0.05, ...
                        'interpolate_time', true, ...
                        'output_rate', 10, ...
                        'boundaries', [0 10], ...
                        'is_circular', false, ...
                        'reaction_params', [0.03 0.055], ...
                        'diffusion_params', [2e-5 1e-5], ...
                        'init_func', {{}}, ...
                        'init_params', [2e-4 1e-4], ...
                        'time_step', 0.5, ...
                        'max_iter', 1e7, ...
                        'tmax', 100, ...
                        'user_data', {{}}, ...
                        'advection_params', [0.02 0.01]);


    % Global structure of a recording/analysis (see ASSET.m)
    case 'mymovie'
      mystruct = struct('correction', [], ...               % Correction between different channels
                        'cortex', [], ...                   % Cortex channel (fluorescence)
                        'data', [], ...                     % Data channel (generic container)
                        'dic', [], ...                      % DIC channel of the experiment & results of the segmentation
                        'eggshell', [], ...                 % Eggshell channel (fluorescence)
                        'experiment', '', ...               % Name of the experiment
                        'markers', [], ...                  % Result of the fluorescent segmentation
                        'metadata', []);

    % Parameters used for the quantification of the signal
    case 'quantification'
      % Retrieve the previously defined structures for the smoothness and data terms
      init = get_struct('data_parameters');
      params = get_struct('smoothness_parameters');
      weights = get_struct('data_parameters');

      params.nhood = 5;
      params.alpha = 0.25;
      params.beta = 0.75;
      params.gamma = 0.75;
      params.delta = 0.5;

      weights.alpha = 0.5;
      weights.beta = 0.75;
      weights.gamma = 0.25;

      init.alpha = 0.45;
      init.beta = 0.25;
      init.gamma = 1e-03;
      init.delta = 0.4;


      mystruct = struct('channel', 'data', ...              % Quantified channel
                        'field', 'cortex', ...              % Quantified field in the previously defined channel
                        'init_params', init, ...           
                        'normalize', 'cytoplasm', ...       % Type of signal normalization
                        'norm_shift', 2, ...                % Distance to quantify the normalization value
                        'params', params, ...
                        'weights', weights, ...
                        'use_ruffles', true, ...            % Quantify along the ruffles ?
                        'resolution', 0.5, ...
                        'pole_threshold', 1/10, ...
                        'kymograph_type', 'projected', ...
                        'window_shape', 'gaussian', ...     % Shape of the quantification window, can either be a filter or a 'fspecial' type
                        'window_params', 0.5, ...           % Parameters required to compute the filter
                        'window_size', 2);                  % Size of (square) the window

    % Parameters of the reference embryo (see carth2normalized.m)
    case 'reference'
      mystruct = struct('axes_length', [25; 15], ...        % Major and minor radii of the ellipse
                        'centers', [0; 0], ...              % Position of the ellipse
                        'orientations', 0, ...              % Tilt (in radians) of the ellipse
                        'index', 0);                        % Identificator for the current ellipse

    % Structure used to store the detected ruffles (see track_ruffles.m)
    case 'ruffles'
      mystruct = struct('bounds', [], ...                   % Bounds of the aperture of the invagination
                        'carth', NaN(1, 2), ...                    % Cartesian position of the invaginations (Nx2)
                        'cluster', [], ...                  % Time cluster containing the invaginations
                        'paths', {{}}, ...                  % Segmentation of the inner part of the invagination
                        'properties', [], ...               % Various properties computed on each invagination
                        'warped', []);                      % Normalized position

    % Parameters used for a segmentation (eggshell and cortex) (see 'segmentations')
    case 'segmentation' 
      % Retrieve the previously defined structures for the smoothness and data terms
      params = get_struct('smoothness_parameters');
      weights = get_struct('data_parameters');

      mystruct = struct('cortex_params', params, ...        % Smoothness for the cortex
                        'cortex_weights', weights, ...      % Data for the cortex
                        'eggshell_params', params, ...      % Smoothness for the eggshell
                        'eggshell_weights', weights, ...    % Data for the eggshell
                        'estimate', [], ...                 % Field to store parameters for the initial elliptical projection
                        'noise', [], ...                    % Field to store parameters to handle noise (filters usually)
                        'safety', 1.2, ...                  % Additional portion projected for safety (see carthesian_coordinate.m)
                        'scoring_func', {{}});              % Function handle for the scoring functions (first:eggshell, second:cortex)

    % Structure containing all the information for all the segmentations (including parameter values)
    % Start by looking into segment_movie.m and dynamic_programming.m
    case 'segmentations'
      % Get the basic structure
      segment = get_struct('segmentation');

      % Construct the global structure to store all the different segmentations
      mystruct = get_struct('mymovie');

      % Set the segmentation to the used channels
      mystruct.dic = segment;                               % Parameters to segment DIC images
      mystruct.data = segment;                              % Parameters to segment data images
      mystruct.markers = segment;                           % Parameters to segment fluorescent images
      mystruct.correction = get_struct('conversion');       % Parameters to convert from DIC to fluorescence

      %%### DIC PARAMETERS ###%%

      % Best values for the eggshell
      % This first section of parameters is common for all the functions
      mystruct.dic.eggshell_params.nhood   = 5;      % Neighborhood size
      mystruct.dic.eggshell_params.alpha   = 0.8873; % Prop. of smoothness VS data
      mystruct.dic.eggshell_params.beta    = 0.1813; % Prop. of path VS intensity
      mystruct.dic.eggshell_params.gamma   = 0.9975; % Prop. of dx VS d2x

      % This part is function-specific
      mystruct.dic.eggshell_weights.alpha  = 0.2490; % Prop. of edges VS outside
      mystruct.dic.eggshell_weights.beta   = 0.0240; % Rescaling factor for the edges
      mystruct.dic.eggshell_weights.eta    = 0.4953; % Position of the eggshell between outside & inside

      % Best values for the cortex
      % Same as for the DIC eggshell
      mystruct.dic.cortex_params.nhood    = 5;
      mystruct.dic.cortex_params.alpha    = 0.7773;
      mystruct.dic.cortex_params.beta     = 0.0769;
      mystruct.dic.cortex_params.gamma    = 0.3467;

      mystruct.dic.cortex_weights.alpha   = 0.2090; % Prop. of edges VS rest
      mystruct.dic.cortex_weights.beta    = 0.5326; % Prop. of outside VS intensity
      mystruct.dic.cortex_weights.gamma   = 0.6979; % Rescaling factor for the edges
      mystruct.dic.cortex_weights.delta   = 0.9854; % Prop. of edges VS gap penalty
      mystruct.dic.cortex_weights.epsilon = 0.8751; % Prop. of cortex VS egg intensity

      mystruct.dic.scoring_func = {@weight_egg, ... % Scoring functions used to
                                   @weight_cortex}; % segment DIC images

      %%### DATA PARAMETERS ###%%

      % Best values for the eggshell
      % Same as for the DIC eggshell
      mystruct.data.eggshell_params.nhood   = 5;
      mystruct.data.eggshell_params.alpha   = 0.9243;
      mystruct.data.eggshell_params.beta    = 0.1458;
      mystruct.data.eggshell_params.gamma   = 0.8584;

      % The actual value of this parameter depends on the pixel size which is not 
      % yet known, consequently it will be computed by "set_pixel_size". 
      % For more information on defining parameters proportional
      % to the size of the pixels: help set_pixel_size.
      mystruct.data.eggshell_weights.filt   = ...    % local "Edge-detection" filter
          ['filt = ones(1,2*ceil(1.5 / pixel_size) + 1);' ...
          'filt(1,1:floor(length(filt) / 2)) = -1;' ...
          'filt / sum(abs(filt));'];
      mystruct.data.eggshell_weights.alpha  = 0.375;% Prop. of intensity VS filter

      % Best values for the cortex
      % Same as for the DIC eggshell
      mystruct.data.cortex_params.nhood  = 9;
      mystruct.data.cortex_params.alpha  = 0.3336;
      mystruct.data.cortex_params.beta   = 0.2941;
      mystruct.data.cortex_params.gamma  = 0.6278;

      mystruct.data.cortex_weights.alpha = 0.5228;   % Prop. of intensity VS outside
      mystruct.data.cortex_weights.beta  = 1e-8;     % Threshold for outside
      mystruct.data.cortex_weights.gamma = 1;        % Intensity target value

      mystruct.data.scoring_func = {@intens_filt,... % Segmentation functions
                                    @intens_sum};

      % General filters used in the segmentation which should be proportional to the 
      % size of the images. For more information on defining parameters proportional
      % to the size of the pixels: help set_pixel_size.
      mystruct.data.shrink = 'strel(''disk'', ceil(2 / pixel_size), 0);';
      mystruct.data.noise  = struct('gaussian', '0.15 / pixel_size;', ...
                                         'median', 'ceil(0.5 / pixel_size) * ones(1,2);');

      %%### MARKERS PARAMETERS ###%%

      % Best values for the eggshell
      % Same as for the DIC eggshell
      mystruct.markers.eggshell_params.nhood   = 5;
      mystruct.markers.eggshell_params.alpha   = 0.9243;
      mystruct.markers.eggshell_params.beta    = 0.1458;
      mystruct.markers.eggshell_params.gamma   = 0.8584;

      % The actual value of this parameter depends on the pixel size which is not 
      % yet known, consequently it will be computed by "set_pixel_size". 
      % For more information on defining parameters proportional
      % to the size of the pixels: help set_pixel_size.
      mystruct.markers.eggshell_weights.filt   = ...    % local "Edge-detection" filter
          ['filt = ones(1,2*ceil(1.5 / pixel_size) + 1);' ...
          'filt(1,1:floor(length(filt) / 2)) = -0.5;' ...
          'filt / sum(abs(filt));'];
      mystruct.markers.eggshell_weights.alpha  = 0.0018;% Prop. of intensity VS filter

      % Best values for the cortex
      % Same as for the DIC eggshell
      mystruct.markers.cortex_params.nhood  = 9;
      mystruct.markers.cortex_params.alpha  = 0.3336;
      mystruct.markers.cortex_params.beta   = 0.2941;
      mystruct.markers.cortex_params.gamma  = 0.6278;

      mystruct.markers.cortex_weights.alpha = 0.5228;   % Prop. of intensity VS outside
      mystruct.markers.cortex_weights.beta  = 1e-8;     % Threshold for outside
      mystruct.markers.cortex_weights.gamma = 1;        % Intensity target value

      mystruct.markers.scoring_func = {@intens_filt,... % Segmentation functions
                                       @intens_sum};

      % General filters used in the segmentation which should be proportional to the 
      % size of the images. For more information on defining parameters proportional
      % to the size of the pixels: help set_pixel_size.
      mystruct.markers.shrink = 'strel(''disk'', ceil(2 / pixel_size), 0);';
      mystruct.markers.noise  = struct('gaussian', '0.15 / pixel_size;', ...
                                         'median', 'ceil(0.5 / pixel_size) * ones(1,2);');

    % Structure used to store the smoothness parameters (see 'segmentations')
    case 'smoothness_parameters'
      mystruct = struct('final', [], ...                % Final position used for backtracking DP
                        'init', [], ...                 % Initial position for DP (see dynamic_programming.m)
                        'nhood', 0, ...                 % Neighborhood explored during dynamic programming (nhood pixels on each side)
                        'prohibit', 'none', ...         % Prohibiting particular moves
                        'spawn_percentile', [], ...     % Score used when spawning a new path (as percentile of the previous step)
                        'alpha', 0, ...                 % Weights of the different smoothness terms
                        'beta', 0, ...                  %  "
                        'gamma', 0, ...                 %  "
                        'delta', 0);                    %  "

    % Structure of MATLAB's splines, useful to represent empty splines (help spline)
    case 'spline'
      mystruct = struct('breaks', [], ...
                        'coefs', [], ...
                        'dim', 0, ...
                        'form', '', ...
                        'order', 0, ...
                        'pieces', 0);

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
      mystruct = struct('fusion_thresh', 2.5, ...         % Minimal distance in um to another spot (estimation) before fusion
                        'frame_displacement', 0.5, ...    % Maximal displacement of a spot (in um) between two frames 
                        'frame_window', 3, ...          % Considered number of frames for the gap closing algorithm (see track_spots.m)
                        'gap_function', @relative_distance, ... % Function used to measure the gap-closing weight
                        'joining_function', @merging_distance, ... % Same but for the joinging weight
                        'splitting_function', @splitting_distance, ... % For the splitting weight
                        'linking_function', @spot_similarity, ... % And for the frame-to-frame linking 
                        'max_size', 1, ...              % Maximal size (in um) of the spots
                        'max_particles', 5, ...
                        'max_iterations', 20, ...
                        'intensity_thresh', 30, ...
                        'fit_intensities', 'separate', ...
                        'iteration_threshold', 1e-10, ...
                        'noise_thresh', 1);             % Threshold used to remove the nosie (see imatrou.m)

    % Structure used to store the parameters required to compute the temperature of the posterior decoding (see find_temperature.m)
    case 'temperatures'
      mystruct = struct('aim_transitions', 1/6, ...     % The temperature of the transitions is based on the mean translation 
                        'aim_emissions', 1/50, ...      % The one for the emission is based on the average width of the posterior path
                        'step_thresh', 1e-5, ...        % Stop criterion in case two consecutive dichotomies have the same value (stalled)
                        'thresh', 1e-3);                % Stop criterion for the dichotomy (diff. to the aims)

    % Structure used to store the timing of the cell cycle
    case 'timing'
      mystruct = struct('cytokinesis', NaN, ...        % Frame of cytokinesis onset
                        'pronuclear_meeting', NaN, ... % Frame of pronuclear meeting
                        'pseudocleavage', NaN);        % Frame of max pseudocleavage

    % Structure containing the data from the manual trackings (see import_trackings.m)
    case 'tracking'
      mystruct = struct('child',[], ...                 % Children structures containing either groups of tracking or tracking files
                        'average', [],  ...
                        'errors',[], ...                % Compute error to the average tracking
                        'expr','', ...                  % Regular expression used to find the children files for import
                        'name','', ...                  % Name of the manual trackings extracted from the filenames
                        'reference', get_struct('reference'));  % Reference embryo used to project the different trackings

    % General structure to store all the manual trackings
    case 'trackings'
      % Get the substructure to store manual trackings
      tracking = get_struct('tracking');

      % Get the structure to store the various channels of a movie
      mystruct = get_struct('mymovie');

      % Set the trackings to the used channels
      mystruct.dic = tracking;
      mystruct.markers = tracking;

    % The structure containing the information to convert into the absolute coordinate system
    case 'warper'
      ref = get_struct('reference');
      orig = ref;
      orig.centers = NaN(2,1);
      orig.axes_length = NaN(2,1);

      mystruct = struct('original', orig,...
                        'reference', ref, ...
                        'warp', []);

    % The structure containing the information to correlate the Z-axis with the others
    case 'z-correlation'
      mystruct = struct('bkg', 18.3835, ...
                        'long_axis', 0.0098, ...
                        'short_axis', -0.2977);

    % If the required type of structure has not been implemented, return an empty one
    otherwise
      mystruct = struct();
  end

  % Repeat the structure to fit the size (nstruct can be multi-dimensional)
  mystruct = repmat(mystruct,nstruct);

  return;
end
