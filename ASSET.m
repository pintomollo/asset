function [mymovie,trackings] = ASSET(varargin)
% ASSET Algorithm for the Segmentation and the Standardization of C. elegans Time-lapse recordings
%   This is the main function for the automated analysis of C. elegans embryos.
%
%   MYMOVIE = ASSET(MYMOVIE, TRACKINGS, OPTS, 'field1', VALUE1, ...) runs an analysis on the
%   recording MYMOVIE using the manual trackings TRACKINGS and the options OPTS. Following
%   arguments are field/VALUE pairs of options (see get_struct.m or get_struct('ASSET')).
%   The returned MYMOVIE is the updated version of the corresponding input. All three
%   structures are automatically saved in a .mat file which name corresponds to the experiment.
%
%   MYMOVIE = ASSET('experiment', ...) loads MYMOVIE, TRACKINGS and OPTS from a previous
%   experiment which results were stored in the MAT-file 'experiment.mat'. Additional options
%   can still be provided for the new analysis.
%
%   MYMOVIE = ASSET runs an analysis by prompting for the recording(s) and the tracking files (if needed).
%   Additional options can provided.
%
%   [MYMOVIE, TRACKINGS] = ASSET(...) also returns the updated TRACKINGS.
%   
%   OPTS most useful fields (default value):
%     - application         ('')/'ruffling'/...     Additional(s) analysis based on the segmentation, can
%                                                   be a cell array (see perform_application.m)
%     - auto_save           true/(false)            Save after every step of the analysis
%     - ccd_pixel_size      (6.45)/...              Camera's pixel size to convert pixels to µm
%     - config_file         ('')/'file.txt'         File containing all the field/VALUE pairs
%                                                   (see Config/example.txt)
%     - magnification       (63)/...                Microscope's magnification to convert pixels to µm
%     - overwrite           (true)/false            Saves the results in the same MAT-file as before
%     - recompute           true/(false)            Recompute the whole analysis
%     - segmentation_type   ('dic')/'markers'/...   Channel(s) used to perform the segmentation, can be
%                                                   a cell array (see segment_movie.m)
%     - verbosity           0/1/(2)/3               Level of output of ASSET. 0 = null, 1 = text only,
%                                                   2 = text/gui, 3 = full (with plots)
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 10.12.2010

  % Check whether ASSET needs to be installed
  install_ASSET

  % Parsing the variable inputs
  [mymovie, trackings, opts] = parse_input(varargin{:});

  % This can occur only if we ran several analysis in a row.
  % So the work was done and we can stop.
  if (isempty(mymovie))
    return;
  end

  % Temporary variables
  expr_name = '';
  is_original_file = true;

  % Debug mode breakpoints
  if (opts.debug)
    display('Inputs parsed and parameters loaded');
    keyboard;

  % Progress display
  elseif (opts.verbosity > 0)
    display('Parameters loaded');
  end

  % We need manual trackings only if we measure the performances of ASSET
  % or if we want to perform machine learning
  if (opts.measure_performances || ~strncmp(opts.do_ml, 'none', 4))
    
    % Import the trackings from the .shapes files
    [trackings, opts, updated] = import_trackings(trackings, opts);

    % Auto-save
    if (updated && opts.auto_save)

      % Use whichever name is available
      tmp_name = '';
      if (~isempty(mymovie.experiment))
        tmp_name = mymovie.experiment;
      elseif (~isempty(trackings.experiment))
        tmp_name = trackings.experiment;
      end

      % Update the experiment name if we cannot overwrite
      if (~opts.overwrite)
        tmp_name = get_new_name([tmp_name '(\d+)\.mat']);

        % We do not keep the extension
        tmp_name = tmp_name(1:end-4);
        
        % Remember we already created a new file
        is_original_file = false;
      end

      % Save if we know where to
      if (~isempty(tmp_name))
        save(tmp_name, 'mymovie', 'trackings','opts');
      end
    end

    % Debug mode breakpoints
    if (opts.debug)
      display('Manual trackings imported');
      keyboard;

    % Progress display
    elseif (opts.verbosity > 0)
      display('Trackings loaded');
    end
  end

  % Import the recordings
  mymovie = open_movie(mymovie, expr_name, opts);

  % Debug mode breakpoints
  if (opts.debug)
    display('Recording imported');
    keyboard;

  % If there is no experiment, something is really wrong, give up
  elseif (isempty(mymovie.experiment))
    error('No recording was properly loaded thus ASSET cannot run');

  % Progress display
  elseif (opts.verbosity > 0)
    display(['Movie ''' mymovie.experiment ''' loaded']);
  end

  % Compute the actual pixel size of the image
  opts = set_pixel_size(opts);

  % Update the experiment name if we cannot overwrite (and if not done before)
  if (~opts.overwrite & is_original_file)
    mymovie.experiment = get_new_name([mymovie.experiment '(\d+)\.mat']);

    % We do not keep the extension
    mymovie.experiment = mymovie.experiment(1:end-4);
  end

  % Runs machine learning algorithms to optimze the segmentation parameters
  if (~strncmp(opts.do_ml, 'none', 4))
    if (opts.verbosity > 0)
      display('Machine learning starts !');
    end

    opts = find_parameters(mymovie, trackings, opts);

    % Debug mode breakpoints
    if (opts.debug)
      display('M-L finished');
      keyboard;

    % Progress display
    elseif (opts.verbosity > 0)
      display('Machien learning has finished');
    end
  end

  % Auto-save
  if (opts.auto_save)
    save(mymovie.experiment, 'mymovie', 'trackings','opts');
  end

  % Segment the movie, the condition is useful if you want to use 'recompute'
  % without re-segmenting the recording
  if (opts.segment)
    [mymovie, updated] = segment_movie(mymovie, opts);

    % Auto-save if something was changed
    if (opts.auto_save && updated)
      save(mymovie.experiment, 'mymovie', 'trackings','opts');
    end

    % Progress display
    if (opts.verbosity > 0)
      display('Movie segmented');
    end
    % Debug mode breakpoints
    if (opts.debug)
      keyboard;
    end
  end

  % Export the movie now that we know more about it
  if (opts.export_movie)
    export_movie(mymovie, opts);

    % Progress display
    if (opts.verbosity > 0)
      display('Movie exported');
    end
    % Debug mode breakpoints
    if (opts.debug)
      keyboard;
    end
  end
  
  % Measuring the performances of ASSET compared to manual tracking
  if (opts.measure_performances)
    [mymovie, trackings] = analyze_segmentation(mymovie, trackings, opts);

    % Auto-save
    if (opts.auto_save)
      save(mymovie.experiment, 'mymovie', 'trackings','opts');
    end

    % Progress display
    if (opts.verbosity > 0)
      display('Precision with respect to ''trackings'' measured');
    % Debug mode breakpoints
    elseif (opts.debug)
      display('Analyze segmentation finished');
      keyboard;
    end
  end

  % Finally perform the application(s)
  if (length(opts.application) > 0)
    mymovie = perform_application(mymovie, opts);

    % Auto-save
    if (opts.auto_save)
      save(mymovie.experiment, 'mymovie', 'trackings','opts');
    end

    % Progress display
    if (opts.verbosity > 0)
      display('Applications computed');
    % Debug mode breakpoints
    elseif (opts.debug)
      display('Applications all performed');
      keyboard;
    end
  end

  % Normalize the movie into the absolute coordinate system
  if (opts.normalize)
    mymovie = carth2normalized(mymovie, opts);

    % Auto-save
    if (opts.auto_save)
      save(mymovie.experiment, 'mymovie', 'trackings','opts');
    end

    % Progress display
    if (opts.verbosity > 0)
      display('Normalization done');
    % Debug mode breakpoints
    elseif (opts.debug)
      display('Data normalized in the absolute coordinate system');
      keyboard;
    end
  end

  % Really save at least once
  if (~opts.auto_save)
    save(mymovie.experiment, 'mymovie', 'trackings','opts');
  end

  % It's over, let's notify it !
  if (opts.verbosity > 1)
    try
      mariosong;
    catch ME
    end
  end

  return;
end

% Sort out the mess due to the variable number of inputs
function [mymovie, trackings, opts] = parse_input(varargin)

  % Initialize the outputs to be sure that we don't return unassigned variables
  mymovie = [];
  trackings = [];
  opts = [];
  
  % Check what we got as inputs
  if (nargin > 0)

    % The first argument should be mymovie, verify using 'experiment'
    if (isfield(varargin{1}, 'experiment'))
      
      % Assign it and delete it from the list
      mymovie = varargin{1};
      varargin(1) = [];

    % Maybe the name of the MAT-file was provided
    elseif (ischar(varargin{1}))
      
      % We allow several files to be analyzed consecutively by using
      % the regular expression metacharacter '*'
      if (~isempty(findstr(varargin{1}, '*')))

        % List all the MAT-files corresponding to the pattern
        datas = dir(varargin{1});

        % Chech whether the file is in a subdirectory
        indx = strfind(varargin{1}, filesep);
        dirpath = '';
        if (~isempty(indx))
          dirpath = varargin{1}(1:indx(end));
        end

        % Run them one after the other one
        for i=1:length(datas)
          ASSET([dirpath datas(i).name], varargin{2:end});
        end

        % Do not forget to stop here !
        return;
      end

      % Now we'll try to load the corresponding MAT-file.
      % So for now on, we suppose that it's really a MAT-file
      % but that we have not loaded it yet.
      done = false;
      real_mat = true;
      waiting_time = 0;

      % We'll try until we loaded it
      while (~done)

        % Accessing a MAT-file which is being accessed by another instance
        % of ASSET produces an error, so we'll catch it.
        try

          % If the name is explicitly a MAT-file which exists, we just load it
          if (strncmp(varargin{1}(end-3:end), '.mat', 4) && exist(varargin{1}, 'file'))
            % Loading will assign the variables
            load(varargin{1});
            % We remove it from the list
            varargin(1) = [];

          % Maybe they forgot the extension ?
          elseif (exist([varargin{1} '.mat'], 'file'))
            load([varargin{1} '.mat']);
            varargin(1) = [];

          % Otherwise it's already a field/VALUE parameter
          else
            real_mat = false;
          end

          % We did our job !
          done = true;

        % Here we intercept the error and will basically wait for
        % some time and try to access it later.
        catch ME

          % If it's a real MAT-file and if we have not already waited
          % too long, then we'll pause again
          if (real_mat && waiting_time < 600)

            % We wait a random amount of time to avoid the (unlikely)
            % case where all the instances wake up at the same time.
            nsecs = rand(1) * 30;
            waiting_time = waiting_time + nsecs;

            % Always display this
            disp(['MAT file is currently used, trying again in ' str2double(nsecs) 's']);
            % Wait
            pause(nsecs);
          else

            % There is a different type of error
            error(ME);
          end
        end
      end
    end
  end

  % If mymovie was not provided, we use the standard structure
  if (isempty(mymovie))
    if (~real_mat & mod(length(varargin), 2) == 1 & ischar(varargin{1}))
      mymovie = varargin{1};
      varargin(1) = [];
    else
      mymovie = get_struct('mymovie');
    end
  end

  % The second possible argument is trackings, it needs to be a
  % structure with a field 'experiment'
  if (~isempty(varargin) && isfield(varargin{1}, 'experiment'))
    % Assign & remove
    trackings = varargin{1};
    varargin(1) = [];
  else

    % Otherwise we assign an empty structure (as it might never be used)
    trackings = get_struct('tracking', 0);
  end

  % Retrieve the latest version of OPTS
  new_opts = get_struct('ASSET');
  new_opts = load_parameters('default_params.txt');

  % Next possible structure has to be opts
  if (~isempty(varargin) && isstruct(varargin{1}))
    opts = varargin{1};
    varargin(1) = [];
  end

  % Correct the structure of opts if need be
  opts = merge_structures(opts, new_opts);

  % Now we check that the parameters were provided in pairs
  npairs = length(varargin) / 2;
  if (npairs ~= floor(npairs))
    error 'Properties pairs must come in PAIRS.';
  end

  % Loop over the pairs of parameters
  for i = 1:npairs

    % If the parameter exists in opts we simply assign it the
    % provided value
    if (isfield(opts, varargin{2*i - 1}))
      opts.(varargin{2*i - 1}) = varargin{2*i};

    % Otherwise we notify its inexistence
    else
      tokens = regexp(varargin{2*i - 1}, '\.', 'split');
      tmp_struct = opts;
      is_valid = true;

      for c = 1:length(tokens)
        if (~isfield(tmp_struct, tokens{c}))
          is_valid = false;
          break;
        else
          tmp_struct = tmp_struct.(tokens{c});
        end
      end

      if (is_valid)
        try
          eval(['opts.' varargin{2*i - 1} ' = varargin{2*i};']);
        catch ME
          is_valid = false;
        end
      end

      if (~is_valid)
        warning(['Property ''' varargin{2*i -1} ''' does not exist. Ignoring']);
      end
    end
  end

  % Try to load the parameters from a configuration file
  opts = load_parameters(opts);

  return;
end
