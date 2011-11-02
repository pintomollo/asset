function [mymovie] = open_movie(mymovie, opts)
% OPEN_MOVIE prepares a recording such that it can be analyzed by ASSET.
%
%   [MYMOVIE] = OPEN_MOVIE(NAME, OPTS) loads the movie NAME into MYMOVIE
%   using the provided options from OPTS. If NAME is empty, the user will be
%   prompted to choose the adequate file.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 19.05.2011

  % Check whether mymovie does need to load the movie.
  % We take into account: 1. mymovie null structure
  %                       2. mymovie contains a null DIC field
  %                       3. either the eggshell or the cortex field is empty
  if(~isstruct(mymovie) | length(mymovie) == 0 | ...
    (isempty(mymovie.dic) & isempty(mymovie.eggshell & isempty(mymovie.cortex))))

    % Initialization
    curdir = '';
    dirpath = '';

    % We got the file to be loaded as input !
    if (ischar(mymovie))

      % Chech whether the file is in a subdirectory
      indx = strfind(mymovie, filesep);

      % If this is the case, separate the path and the name
      if (~isempty(indx))
        dirpath = mymovie(1:indx(end));
        fname = mymovie(indx(end)+1:end);

        % Check whether this is an absolute path or not
        if (dirpath(1) ~= filesep && isempty(strfind(mymovie, ':')))
          dirpath = ['.' filesep dirpath];
        end
      else
        fname = mymovie;
      end
    else
  
      % In case a subfolder name Movies exists, move into it for prompting
      curdir = '';
      if(exist('Movies', 'dir'))
        curdir = cd;
        cd('Movies');
      elseif(exist(['..' filesep 'Movies'], 'dir'))
        curdir = cd;
        cd(['..' filesep 'Movies']);
      end

      % Fancy output
      disp('[Select a movie file]');

      % Prompting the user for the movie file
      [fname, dirpath] = uigetfile({'*.*'}, ['Load a movie']);
    end

    % Return back to our original folder
    if(~isempty(curdir))
      cd(curdir);
    end

    % If no file was selected, stop here
    if (length(fname) == 0  ||  isequal(dirpath, 0))
      disp(['No movie selected']);
      return;
    end  

    % Parse the file name using the provided pattern to identify the other channels
    [tokens, junk] = regexp(fname, opts.file_regexpr, 'tokens');

    name = tokens{1}{1};
    suffix = tokens{1}{2};
    ext = tokens{1}{3};

    % If the user chose a MAT file, load it and stop here
    if (strncmp(ext, '.mat', 4))
      load(fname);
      return;
    end

    % Check if we have the "full" pattern
    if (~isempty(name))

      % List the available files that match the pattern
      files = dir([dirpath name '*' ext]);

      % Parse the different available files, removing the ones that do no fit
      % the pattern, as the introduction of the "*" might give some false hits
      for i = length(files):-1:1
        [tokens, junk] = regexp(files(i).name, opts.file_regexpr, 'tokens');

        if (~strcmp(name, tokens{1}{1}))
          files(i) = [];
        end
      end

    % Otherwise, use directly the name
    else
      name = suffix;
      suffix = '';

      files = struct('name', [name ext]);
    end

    % Count the number of channels and prepare the structure accordingly
    nchannels = length(files);
    channels = get_struct('channel',[1 nchannels]);

    % Retrieve the base directory of the recordings
    base_dir = '';

    % Open each channel separatly. The policy is used to remember the user's choice 
    policy = 0;
    for i=1:nchannels
      channels(i).file = [dirpath files(i).name];
      channels(i).file = relativepath(channels(i).file);

      % Get the directory in which the reocrdings are stored
      if (isempty(base_dir))
        [base_dir] = fileparts(channels(i).file);
      end

      % We convert the provided type into a more handy one
      [channels(i).fname, policy] = convert_movie(channels(i).file, policy, opts);
    end

    % Ask the user to identify the different channels
    mymovie = identify_channels(channels, opts);

    % Use the movie name as the name for the entire experiment
    mymovie.experiment = name;

    % Load the metadata if they exist
    mymovie = parse_metadata(mymovie, base_dir, opts);

    % Rescale the channels (including some filtering)
    mymovie = rescale_movie(mymovie, opts);
  end

  return;
end
