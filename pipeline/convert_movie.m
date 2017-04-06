function [converted_file, metadata, opts] = convert_movie(name)
% CONVERT_MOVIE converts a recording such that it can be tracked properly.
%
%   [CONVERTED] = CONVERT_MOVIE(NAME) converts the movie NAME into CONVERTED
%   which is a single OME-TIFF stack file for convenience. The LOCI command toolbox
%   is used for this process. This allows to work with a single file type later on.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 04.07.2014

  % Initialization
  converted_file = '';
  metadata = '';
  opts = '';
  curdir = '';
  dirpath = '';
  if (nargin == 0)
    name = '';
  end

  % We got the file to be loaded as input !
  if (ischar(name) && ~isempty(name))

    % Chech whether the file is in a subdirectory
    indx = strfind(name, filesep);

    % If this is the case, separate the path and the name
    if (~isempty(indx))
      dirpath = name(1:indx(end));
      fname = name(indx(end)+1:end);

      % Check whether this is an absolute path or not
      if (dirpath(1) ~= filesep && isempty(strfind(name, ':')))
        dirpath = ['.' filesep dirpath];
      end

      % Reconstruct the full path
      fname = fullfile(dirpath, fname);
    else
      fname = name;
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
    fname = fullfile(dirpath, fname);
  end

  % Return back to our original folder
  if(~isempty(curdir))
    cd(curdir);
  end

  % If no file was selected, stop here
  if (isempty(fname)  ||  isequal(dirpath, 0))
    disp(['No movie selected']);
    return;
  end

  % Check the type of the movie chosen
  [file_path, filename, ext] = fileparts(fname);

  % If the user chose a MAT file, load it and stop here
  if (strncmp(ext, '.mat', 4))
    load(fname);
    return;
  end

  % This can take a while, so inform the user
  hInfo = warndlg('Parsing metadata, please wait.', 'Preprocessing movie...');

  % Try to identify better metadata
  metadata = find_metadata(fname);

  % Store the resulting metadata
  [metadata, opts] = parse_metadata(metadata);

  keyboard

  % Close the info
  if (ishandle(hInfo))
    delete(hInfo);
  end

  % Did we get a valid list of metadata files ?
  if (~isempty(metadata.files) && exist(fullfile(file_path, metadata.files{1}), 'file'))

      %%% NEED TO AVOID CONVERTING TIFF FILES
      %
      % Easy way of telling if we have a TIFF file already
      %%[file_path, filename, ext] = fileparts(myrecording);
      %%if (~strncmp(ext, '.tif', 4) && ~strncmp(ext, '.tiff', 5))
      %%  myrecording = convert_movie(myrecording, false);
      %%end


    % This also takes quite some time, so warn the user
    hInfo = warndlg('Converting to TIFF, please wait...', 'Converting movie...');

    [junk, dirname, junk] = fileparts(file_path);
    converted_file = cell(length(metadata.channels), 1);

    % We initially do not know what to do
    answer = 0;

    % Get ready to combine the various files together
    metadata.files = fullfile(file_path, metadata.files);
    for i = 1:length(metadata.channels)
      files = metadata.files(i,:,:);
      files = files(:).';

      % Create a new name
      tmp_name = [dirname '_' metadata.channels{i} '.tiff'];
      converted_file{i} = fullfile(file_path, strrep(tmp_name, ' ', '-'));

      % If the file already exists, we ask what to do
      if(exist(converted_file{i},'file'))
        % Creat the fancy name for display (otherwise it thinks they are LaTeX commands)
        [junk, tmp_name, junk] = fileparts(converted_file{i});
        printname = strrep(tmp_name,'_','\_');

        % We do not accept "empty" answers
        while (answer == 0)
          answer = menu(['The TIFF version of ' printname ' already exists, overwrite it ?'],'Yes','No');
        end

        % Act accorindly to the policy
        if (answer == 1)
          delete(converted_file{i});
        end
      end

      % Convert the file
      save_data(converted_file{i}, files, 'uint16');
    end

    % Close the info
    if (ishandle(hInfo))
      delete(hInfo);
    end

  % We do not have valid data
  else
    disp(['No valid movie found']);
    metadata = '';
    opts = '';
  end

  return;
end

function metadata = find_metadata(filename)
% This function tries to identify more suitable metadata. For now
% on, the following metadata are supported:
%   - Leica Application Suite ".las"
%   - Leica Application Suite ".cal.xml"
%   - uManager "metadata.txt"
%   - files manually placed in the "Metadata" folder and named
%     as the recording

  % Initialize the metadata
  metadata = '';

  % Get the folder in which the file is contained
  [file_path, file_name, file_ext] = fileparts(filename);

  % Check if the file exists
  if (exist(fullfile(file_path, '.las'), 'file'))

    % Load it !
    metadata = fileread(fullfile(file_path, '.las'));

  % The other LAS
  elseif (exist(fullfile(file_path, '.Metadata'), 'dir'))
    if (exist(fullfile(file_path, '.Metadata', [file_name file_ext '.cal.xml']), 'file'))
      metadata = fileread(fullfile(file_path, '.Metadata', [file_name file_ext '.cal.xml']));
    else
      indx = strfind(file_name, '_');
      if (~isempty(indx))
        file = regexpdir(fullfile(file_path, '.Metadata'), [file_name(1:indx(end)-1) '(\..*)?\.cal\.xml' ]);
        if (length(file)==1)
          metadata = fileread(file{1});
        end
      end
    end

  % For uManager
  elseif (exist(fullfile(file_path, 'metadata.txt'), 'file'))
    metadata = fileread(fullfile(file_path, 'metadata.txt'));

  % For manually edited files
  elseif (exist(fullfile(pwd, 'Metadata', [filename '.txt']), 'file'))
    metadata = fileread(fullfile(pwd, 'Metadata', [filename '.txt']));

  % As a last resort, we'll try to use the LOCI tool if installed
  else

    % Solution to Java heap space errors:
    % https://www.mathworks.com/matlabcentral/answers/92813

    % We need the absolute path for Java to work properly
    filename = absolutepath(filename);

    if (ispc)
      cmd_name = ['"' filename '"'];
    else
      cmd_name = strrep(filename,' ','\ ');
    end

    % Look for the LOCI command line tool
    curdir = pwd;
    cmd_path = which('bfconvert.bat');

    % If it's installed, run it
    if (~isempty(cmd_path))
      [mypath, junk] = fileparts(cmd_path);

      % Move to the correct folder
      cd(mypath);

      % And call the LOCI utility to extract the metadata
      if (ispc)
        [res, metadata] = system(['showinf.bat -stitch -nopix -nometa -omexml-only ' cmd_name]);
      else
        [res, metadata] = system(['./showinf -stitch -nopix -nometa -omexml-only ' cmd_name]);
      end

      % Check if an error occured, remove the metadata
      if (res ~= 0)

        % A issue with the versions of Java...
        if (~isempty(strfind(metadata, 'JNI')))
          warning('ASSET:javaVersion', 'MATLAB appears to be running on a different version of java than the one installed on your computer, which prevents LOCI from running.\nTo fix this issue, see https://www.mathworks.com/matlabcentral/answers/103056')
        end

        metadata = '';

      % Remove the unused lines of data
      else
        metadata = regexprep(metadata, '^[^=:]*\n', '');
      end

      cd(curdir);
    end
  end

  return;
end


