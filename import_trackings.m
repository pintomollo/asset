function [trackings, opts, updated] = import_trackings(trackings, opts)
% IMPORT_TRACKING loads the content of manual tracking files (.shapes) into the
%   corresponding data structure (get_struct('trackings')). It also accepts MAT-files
%   in which 'trackings' structures have been stored.
%
%   [TRACKINGS, OPTS] = IMPORT_TRACKINGS(TRACKINGS, OPTS) imports the files listed in
%   OPTS.TRACKINGS into TRACKINGS. If no file is specified, it prompts the user for it.
%
%   [...] = IMPORT_TRACKINGS() uses the default values for both OPTS and TRACKINGS.
%   If only one is provided, the default value will be used for the missing argument.
%
%   [TRACKINGS, OPTS, EXPERIMENT, UPDATED] = IMPORT_TRACKINGS(...) returns the experiment
%   name EXPERIMENT extracted from the file names, as well as the UPDATED flag which is true
%   when TRACKINGS was modified.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 15.12.2010

  % Some basic input checking to handle missing arguments
  if (nargin == 1)
    % We identify trackings through this field
    if (isfield(trackings, 'experiment'))
      opts = get_struct('ASSET');
    elseif (isstruct(trackings))
      opts = trackings;
      trackings = [];
    else
      trackings = [];
      opts = get_struct('ASSET');
    end
  elseif (nargin == 0)
    trackings = [];
    opts = get_struct('ASSET');
  end
  was_empty = false;

  % Default value for the return values
  updated = false;
  experiment = '';

  % We extract the file names
  files = opts.trackings;
  if (~iscell(files))
    files = {files};
  end

  % 'all' is a keyword which means both DIC and Markers
  if (strncmp(opts.segmentation_type, 'all', 3))
    types = {'dic', 'markers'};
  else
    types = {opts.segmentation_type};
  end

  % Get the required number of segmentation types
  ntypes = length(types);

  % Get the basic structure for every type
  new_trackings = get_struct('trackings');
  tracking = get_struct('tracking');

  % We check whether trackings is already up-to-date.
  % Obviously if it's empty we have some work to do.
  if (isempty(trackings))
    trackings = new_trackings;
    was_empty = true;

  % Otherwise, we simply verify whether its fields corresponding to the 
  % segmentation types were assigned previously.
  elseif (isstruct(trackings))
    
    for t = ntypes
      % If the field do not exist or was not assigned, we have some work to do !
      if (~isfield(trackings, types{t}))
        trackings.(types{t}) = tracking;
      end
    end
  end

  % Required for finding missing frames of manual tracking
  max_frames = 0;

  % We check the dimensionality of the provided files
  [file_types, file_groups, file_ids] = size(files);
  if (file_types > ntypes)
    warning('More tracking files than segmentation types were provided. Additional files will be ignored');
  end

  % We loop over the files to retrieve the data
  for t = 1:file_types
    tmp_type = types{t};

    % We create a new regular expression representing each group of files and prepare 
    % the tracking structure
    if (~isfield(new_trackings, tmp_type))
      new_trackings.(types{t}) = tracking;
    end
    new_trackings.(tmp_type).expr = '';
    new_trackings.(tmp_type).child = get_struct('tracking', [1 file_groups]);

    for g = 1:file_groups

      new_trackings.(tmp_type).child(g).child = get_struct('file', [1 file_ids]);

      for i = 1:file_ids

        % If the name is empty, we cannot do anything
        if (isempty(files{t, g, i}))
          break;
        end

        % Extract the "name", "suffix" and extension of the filename
        [tokens, junk] = regexp(files{t, g, i}, opts.file_regexpr, 'tokens');

        name = tokens{1}{1};
        suffix = tokens{1}{2};
        ext = tokens{1}{3};

        % If a MAT-file is specify, we try to extract the corresponding files from the tracking structure
        if (~isempty(ext) & strncmp(ext, '.mat', 4))
          tmp_trackings = load(files{t, g, i});

          % We merge both structure to avoid loosing any information
          new_trackings = merge_structures(new_trackings, tmp_trackings.trackings, {'expr', 'fname'});

        % Otherwise it's a plain file
        else

          % Here we try to identify the regular expression
          if (isempty(new_trackings.(tmp_type).child(g).expr))
            
            % The file might have no suffix
            if (isempty(name))
              expr = [suffix ext];
              name = expr;
            else
              expr = [name '*' ext];
            end

            % We copy the useful information
            new_trackings.(types{t}).child(g).name = name;
            new_trackings.(types{t}).child(g).expr = expr;
          end

          % And we copy the filename
          new_trackings.(types{t}).child(g).child(i).fname = files{t}{g}{i};
        end
      end
    end
  end

  % Now we merge the new with the old trackings
  [trackings, is_same] = merge_structures(trackings, new_trackings, {'expr', 'fname'});

  % If nothing changed and we do not need to recompute everything, exit
  if (~was_empty & is_same & ~opts.recompute)
    return;
  end

  % Now we prompt and sort the resulting structure
  for t = 1:length(types)
    % First get the name of the experiment based on the filenames
    experiment = get_name(trackings.(types{t}), '');

    % Then we sort the different groups by name
    child_names = {};
    for i=1:length(trackings.(types{t}).child)
      child_names = [child_names trackings.(types{t}).child(i).name];
    end

    % We we have all the names, we can just reorder them by permutation
    [junk, indxs] = sort(child_names);
    trackings.(types{t}).child = trackings.(types{t}).child(indxs);

    % If verbosity is 'big' enough, we can prompt
    if (opts.verbosity > 1)

      disp(['[Select your ' experiment ' ' types{t} ' tracking files]'])

      % Display the input dialog
      trackings.(types{t}) = input_trackings(trackings.(types{t}), [types{t} ' (' experiment ')'], opts);

      % Update the experiment name if needed
      experiment = get_name(trackings.(types{t}), experiment);
    end

    % Now we finally load the trackings
    [trackings.(types{t}), nframes] = load_trackings(trackings.(types{t}), opts);

    % We keep track of the total number of frames to check later if some are missing
    if (nframes > max_frames)
      max_frames = nframes;
    end
  end

  % We check whether some tracking files are missing in some groups and we
  % buffer them accordingly so that every group has the same number of frames
  trackings = check_trackings(trackings, max_frames, opts.verbosity);

  % Save the name and exit
  trackings.experiment = experiment;
  updated = true;

  return;
end

function name = get_name(mystruct, name)
% This function extracts the common substring from various names in order
% to deduce the 'experiment' name

  % We loop over all the childrens to browse their names
  for i = 1:length(mystruct.child)

    % If they have a 'name' field, we try to find the common part with the 
    % provided name string
    if (isfield(mystruct.child(i), 'name') & ~isempty(mystruct.child(i).name))
      if (isempty(name))
        name = mystruct.child(i).name;
      else
        [name, end_name] = common_substring(name, mystruct.child(i).name);

        % Use the end of the name instead
        if (isempty(name))
          name = end_name;
        end
      end
    end
  end

  % Remove heading and trailing non-alphanumeric symbols which might cause
  % problems when saving the .mat file
  name = regexprep(name, '^\W+|\W+$', '');

  return;
end

function trackings = check_trackings(trackings, max_frames, verbosity)
% This function parses a tracking structure to make sure that no frame is missing.
% Moreover it will add empty 'frames' where they are missing to avoid indexing problems.

  % We need to check at which level of the structure we are
  if (isfield(trackings, 'child'))

    % Here we found some children structures so we recursively parse them as they
    % should contain the filenames at some level
    nchild = length(trackings.child);
    for i=1:nchild
      trackings.child(i) = check_trackings(trackings.child(i), max_frames, verbosity);
    end

  % Here we have reached a file structure and thus can check for the missing frames
  elseif (isfield(trackings, 'shapes'))

    [ngroups, nframes] = size(trackings.shapes);

    % In case of missing frames at the end of the tracking, we would miss them and could have
    % indexing problems so we substitute them with empty cells
    if (nframes < max_frames)
      trackings.shapes = [trackings.shapes cell([ngroups, max_frames - nframes])];
      nframes = max_frames;
    end
    
    % We loop over every frame and announce the missing ones
    for g = 1:ngroups
      for i = 1:nframes
        if (isempty(trackings.shapes{g, i}) & verbosity > 0)
          warning(['File ' trackings.fname ', tracking of ' trackings.groups{g} ' in frame ' num2str(i) ' is missing']);
        end
      end
    end

  % Otherwise, we try to look recursively in every structure field
  else

    % Get the existing fields
    fields = fieldnames(trackings);

    % And loop over the structure ones
    for f = 1:length(fields)
      if (isstruct(trackings.(fields{f})))
        trackings.(fields{f}) = check_trackings(trackings.(fields{f}), max_frames, verbosity);
      end
    end
  end

  return;
end
