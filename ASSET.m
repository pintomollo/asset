function [myrecording, opts] = ASSET(varargin)
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
%     - ccd_pixel_size      (6.45)/...              Camera's pixel size to convert pixels to �m
%     - config_file         ('')/'file.txt'         File containing all the field/VALUE pairs
%                                                   (see Config/example.txt)
%     - magnification       (63)/...                Microscope's magnification to convert pixels to �m
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

  % General error handling
  try
  working = true;

  % Check whether ASSET needs to be installed
  install_ASSET

  % Parsing the variable inputs
  [myrecording, opts] = parse_input(varargin{:});

  % This can occur only if we ran several analysis in a row.
  % So the work was done and we can stop.
  if (isempty(myrecording))
    return;
  end

  % Temporary variables
  is_original_file = true;

  % Debug mode breakpoints
  if (opts.debug)
    display('Inputs parsed and parameters loaded');
    keyboard;

  % Progress display
  elseif (opts.verbosity > 0)
    display('Parameters loaded');
  end

  % Check if we do need to load a movie
  if (~isstruct(myrecording))
    if (ischar(myrecording))

      % Easy way of telling if we have a TIFF file already
      [file_path, filename, ext] = fileparts(myrecording);
      if (~strncmp(ext, '.tif', 4) && ~strncmp(ext, '.tiff', 5))
        myrecording = convert_movie(myrecording, false);
      end
      [myrecording, opts] = inspect_recording(myrecording, opts, batch_mode);
    else
      [myrecording, opts] = inspect_recording();
    end

    % Debug mode breakpoints
    if (opts.debug)
      display('Recording imported');
      keyboard;

    % If there is no experiment, something is really wrong, give up
    elseif (~isstruct(myrecording) || ~isfield(myrecording, 'experiment') || isempty(myrecording.experiment))
      error('No recording was properly loaded thus ASSET cannot run');

    % Progress display
    elseif (opts.verbosity > 0)
      display(['Movie ''' myrecording.experiment ''' loaded']);
    end
  end

  % Update the experiment name if we cannot overwrite (and if not done before)
  if (~opts.overwrite & is_original_file)
    myrecording.experiment = get_new_name([myrecording.experiment '(\d+)\.mat']);

    % We do not keep the extension
    myrecording.experiment = myrecording.experiment(1:end-4);
  end

  % Filter the data channels if need be, still save in just in case it is modified
  save([myrecording.experiment '.mat'], 'myrecording', 'opts');
  [myrecording, opts] = filter_channels(myrecording, opts);

  % Compute the actual pixel size of the image
  opts = set_pixel_size(opts);
  save([myrecording.experiment '.mat'], 'myrecording', 'opts');

  % Segment the movie, the condition is useful if you want to use 'recompute'
  % without re-segmenting the recording
  if (opts.segment)
    [myrecording, updated] = segment_movie(myrecording, opts);
    save([myrecording.experiment '.mat'], 'myrecording', 'opts');

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
    export_movie(myrecording, opts);

    % Progress display
    if (opts.verbosity > 0)
      display('Movie exported');
    end
    % Debug mode breakpoints
    if (opts.debug)
      keyboard;
    end
  end

  % Finally perform the application(s)
  if (length(opts.application) > 0)
    myrecording = perform_application(myrecording, opts);
    save([myrecording.experiment '.mat'], 'myrecording', 'opts');

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
    myrecording = carth2normalized(myrecording, opts);
    save([myrecording.experiment '.mat'], 'myrecording', 'opts');

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
  opts.recompute = false;
  save([myrecording.experiment '.mat'], 'myrecording', 'opts');

  % Catch the error overall
  catch
    err = lasterror();
    if (exist('myrecording', 'var') & isfield(myrecording, 'experiment') & ~isempty(myrecording.experiment))
      warning(['Error during the analysis of ' myrecording.experiment ':']);
    else
      warning(['Error during the analysis:']);
    end
    print_all(err);
    working = false;

    close all hidden force;
  end

  % It's over, let's notify it !
  if (~exist('opts', 'var') | opts.verbosity > 1)
    try
      mariosong(working);
    catch
    end
  end

  return;
end

% Sort out the mess due to the variable number of inputs
function [myrecording, opts] = parse_input(varargin)

  % Initialize the outputs to be sure that we don't return unassigned variables
  myrecording = [];
  opts = [];
  real_mat = true;
  reset_opts = false;

  % Check what we got as inputs
  if (nargin > 0)

    % The first argument should be myrecording, verify using 'experiment'
    if (isfield(varargin{1}, 'experiment'))

      % Assign it and delete it from the list
      myrecording = varargin{1};
      varargin(1) = [];

    % Maybe several names were provided
    elseif (iscell(varargin{1}) && (mod(length(varargin), 2)==1 || isstruct(varargin{2})))

      files = varargin{1};
      for i=1:length(files)
        ASSET(files{i}, varargin{2:end});
      end

      % Do not forget to stop here !
      return;

    % Maybe the name of the MAT-file was provided
    elseif (ischar(varargin{1}) && (mod(length(varargin), 2)==1 || isstruct(varargin{2})))

      fname = varargin{1};
      fname = absolutepath(fname);
      [fpath, name, ext] = fileparts(fname);
      files = regexpdir(fpath, [name ext], false);

      if (numel(files) == 0)
        files = regexpdir(pwd, [name ext '.mat'], false);
      end

      if (numel(files) > 1)
        for i=1:length(files)
          ASSET(files{i}, varargin{2:end});
        end

        % Do not forget to stop here !
        return;
      elseif (numel(files) == 1)
        varargin{1} = files{1};
      elseif (~exist(varargin{1}, 'file'))
        warning(['No file could be identified with the pattern: ' varargin{1}]);
        return;
      end

      % Now we'll try to load the corresponding MAT-file.
      % So for now on, we suppose that it's really a MAT-file
      % but that we have not loaded it yet.
      done = false;
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
            % Try closing the possible hanging files
            fclose('all');
          else

            % There is a different type of error
            rethrow(ME);
          end
        end
      end
    end
  end

  % If myrecording was not provided, we use the standard structure
  if (isempty(myrecording))
    if (~real_mat && mod(length(varargin), 2) == 1 && ischar(varargin{1}))
      myrecording = varargin{1};
      varargin(1) = [];
    else
      myrecording = get_struct('myrecording');
    end
  end

  % Retrieve the latest version of OPTS
  new_opts = load_parameters('default_params.txt');

  % Next possible structure has to be opts
  if (~isempty(varargin) && isstruct(varargin{1}))
    opts = varargin{1};
    varargin(1) = [];
  end

  % Correct the different structures of opts if need be
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

    % Otherwise we might be accessing a subfield
    else
      % Find the dots indicating a subfield
      tokens = regexp(varargin{2*i - 1}, '\.', 'split');
      tmp_struct = opts;
      is_valid = true;

      % Check if the whole structure exists
      for c = 1:length(tokens)
        if (~isfield(tmp_struct, tokens{c}))
          is_valid = false;
          break;
        else
          tmp_struct = tmp_struct.(tokens{c});
        end
      end

      % Assign the value
      if (is_valid)
        try
          eval(['opts.' varargin{2*i - 1} ' = varargin{2*i};']);
        catch ME
          is_valid = false;
        end
      end

      % Otherwise we notify its inexistence
      if (~is_valid)
        if (strncmp(varargin{2*i-1}, 'reset', 5))
          reset_opts = varargin{2*i};
        else
          warning(['Property ''' varargin{2*i -1} ''' does not exist. Ignoring']);
        end
      end
    end
  end

  if (reset_opts)
    opts = new_opts;
  end

  % Try to load the parameters from a configuration file
  opts = load_parameters(opts);

  return;
end
