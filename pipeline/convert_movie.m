function [converted_file, metadata, opts] = convert_movie(name, force_do_merge)
% CONVERT_MOVIE converts a recording such that it can be tracked properly.
%
%   [CONVERTED] = CONVERT_MOVIE(NAME) converts the movie NAME into CONVERTED
%   which is a single OME-TIFF stack file for convenience. The LOCI command toolbox
%   is used for this process. This allows to work with a single file type later on.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 04.07.2014

      %%% NEED TO AVOID CONVERTING TIFF FILES
      %
      % Easy way of telling if we have a TIFF file already
      %%[file_path, filename, ext] = fileparts(myrecording);
      %%if (~strncmp(ext, '.tif', 4) && ~strncmp(ext, '.tiff', 5))
      %%  myrecording = convert_movie(myrecording, false);
      %%end


  % Initialization
  converted_file = '';
  metadata = '';
  opts = '';
  curdir = '';
  dirpath = '';
  do_merge = false;
  if (nargin == 0)
    name = '';
    force_do_merge = false;
  elseif (nargin == 1)
    if (islogical(name))
      do_merge = name;
      force_do_merge = true;
      name = '';
    else
      force_do_merge = false;
    end
  elseif (nargin == 2)
    do_merge = force_do_merge;
    force_do_merge = true;
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
  metadata = find_metadata(fname, get_struct('metadata'));

  % Store the resulting metadata
  [metadata, opts] = parse_metadata(metadata);

  % Close the info
  if (ishandle(hInfo))
    delete(hInfo);
  end

  % Did we get a valid a valid list of metadata files ?
  if (~isempty(metadata.files) && exist(fullfile(file_path, metadata.files{1}), 'file'))

    % This also takes quite some time, so warn the user
    hInfo = warndlg('Converting to OME-TIFF, please wait...', 'Converting movie...');

    [junk, dirname, junk] = fileparts(file_path);
    converted_file = cell(length(metadata.channels), 1);

    % We initially do not know what to do
    answer = 0;
    if (do_merge)
      answer = 2;
    end

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
          save_data(converted_file{i}, files, 'uint16');
        end
      else
        save_data(converted_file{i}, files, 'uint16');
      end
    end

    % Close the info
    if (ishandle(hInfo))
      delete(hInfo);
    end
  else

    % Traditionally, channels of an experiment were named the same, but for
    % the last "_channel", so try to identify them
    hits = strfind(filename, '_');
    if (~isempty(hits))
      others = dir(fullfile(file_path, [filename(1:hits(end)) '*']));

      % If may have found more files
      if (length(others) > 1)
        new_files = {};
        for i=1:length(others)
          new_files(end+1) = {others(i).name};
        end

        % Ask a confirmation from the user
        answer = questdlg([{'The following files appear to belong to the same experiment:'}; ...
                            {''}; new_files(:); {''}; {'Open all of them ?'}], ...
                            'Open multiple files');

        % If he agrees, load all files
        if (strcmp(answer, 'Yes'))
          fname = new_files;
        end
      end
    end

    % We convert the provided type into a more handy one using LOCI
    if (iscell(fname))
      converted_file = fname;
      for i=1:length(fname)
        converted_file{i} = bftools_convert(fullfile(file_path, fname{i}), force_do_merge, do_merge);
      end
    else
      converted_file = bftools_convert(fname, force_do_merge, do_merge);
    end
  end

  return;
end

% Here we'll call the LOCI Bio-formats toolbox to convert anything into our
% favorite OME-TIFF format
function [newfile] = bftools_convert(fname, forced, do_merge)

  % We need the absolute path for Java to work properly
  fname = absolutepath(fname);
  if (exist(fname, 'file') ~= 2)
    error('Tracking:BadFile', ['File ' fname ' does not exist']);
  end

  if (ispc)
    cmd_name = ['"' fname '"'];
  else
    cmd_name = strrep(fname,' ','\ ');
  end

  % Split the filename
  [file_path, filename, ext] = fileparts(fname);

  % If we have an AVI recording try using the build-in reader, or FFMPEG to extract the frames
  use_tmp_folder = false;
  if (strncmpi(ext, '.avi', 4))
    % Try using the built-in video reader
    newfname = movie2tiff(fname, forced);

    % Nothing changed, something went wrong
    if (length(newfname) == length(fname))
      if (ispref('ffmpeg', 'exepath'))

        % We need to store the original names for later
        orig_name = filename;
        orig_path = file_path;
        orig_cmd_name = cmd_name;

        % We will extract everything in a subdirectory
        use_tmp_folder = true;
        tmp_folder = absolutepath(get_new_name('tmp_folder(\d+)', file_path));
        mkdir(tmp_folder);

        % Extract the frames as JPEG
        new_fname = fullfile(tmp_folder, 'tmp_img%d.jpg');
        if (ispc)
          cmd_name_new = ['"' new_fname '"'];
        else
          cmd_name_new = strrep(new_fname,' ','\ ');
        end

        % This can take a while, so inform the user
        hInfo = warndlg('Converting AVI using FFMPEG, please wait...', 'Converting movie...');

        % Get FFMPEG to actually extract the frames
        [res, info] = system([getpref('ffmpeg', 'exepath') ' -i ' cmd_name ' -y -vf select="eq(pict_type\,PICT_TYPE_I)" -vsync 2 -qscale:v 2 -f image2 ' cmd_name_new]);

        % Close the info
        if (ishandle(hInfo))
          delete(hInfo);
        end

        % Show the error
        if (res~=0)
          error('Tracking:FFMPEG', info);
        end

        % Create the new filename pointing to the frames
        fname = fullfile(tmp_folder, 'tmp_img1.jpg');

      % Advise the installation of FFMPEG
      else
        warndlg({'Converting AVI files works best using FFMPEG.', ...
          'Consider installing this library if the current conversion does not work.', ...
          'Follow the instructions from install_leachi_flow.m to do so.'}, 'Converting movie...');
      end
    else
      % Use the new file instead
      fname = newfname;
    end

    % Update the names
    [file_path, filename, ext] = fileparts(fname);
    if (ispc)
      cmd_name = ['"' fname '"'];
    else
      cmd_name = strrep(fname,' ','\ ');
    end
  end

  % Remove the extension
  name = fullfile(file_path, filename);

  % Look for the LOCI command line tool
  curdir = pwd;
  cmd_path = which('bfconvert.bat');
  if (isempty(cmd_path))
    error('Tracking:lociMissing', 'The LOCI command line tools are not present !\nPlease follow the instructions provided by install_cell_tracking');
  end
  [mypath, junk] = fileparts(cmd_path);

  % This can take a while, so inform the user
  hInfo = warndlg('Parsing metadata, please wait...', 'Converting movie...');

  % Move to the correct folder
  cd(mypath);

  % And call the LOCI utility to extract the metadata
  if (ispc)
    [res, metadata] = system(['showinf.bat -stitch -nopix -nometa ' cmd_name]);
  else
    [res, metadata] = system(['./showinf -stitch -nopix -nometa ' cmd_name]);
  end

  % Delete the information if need be
  if (ishandle(hInfo))
    delete(hInfo);
  end

  % Check if an error occured
  if (res ~= 0)
    cd(curdir);
    error(metadata);
  end

  % Extract the three important informations from the extracted metadata
  format = regexp(metadata, 'file format\s*\[([ -\w]+)\]', 'tokens');
  is_rgb = regexp(metadata, 'RGB\s*=\s*(\w+)', 'tokens');
  file_pattern = regexp(metadata, 'File pattern = ([^\n]*)\n', 'tokens');
  is_grayscale = regexp(metadata, 'SizeC\s*=\s*\d+\s*(\(effectively 1\))', 'tokens');

  % In case of multiple files, regroup them into one single file
  if (forced)
    if (do_merge)
      merge_cmd = '-stitch ';
      if (~isempty(file_pattern) && ~use_tmp_folder)
        orig_pattern = file_pattern{1};
        file_pattern = regexprep(orig_pattern, '<\d+-\d+>', '');
        if (length(orig_pattern{1}) == length(file_pattern{1}))
          file_pattern = '';
        end
      end
    else
      merge_cmd = '';
      file_pattern = '';
    end
  else
    merge_cmd = '-stitch ';
    do_merge = use_tmp_folder;
    if (~isempty(file_pattern) && ~use_tmp_folder)
      orig_pattern = file_pattern{1};
      file_pattern = regexprep(orig_pattern, '<\d+-\d+>', '');

      % In case we did not delete the pattern, it means there was nothing to merge !
      if (length(orig_pattern{1}) == length(file_pattern{1}))
        file_pattern = '';
      else

        % Otherwise, make sure they should be merged !!
        answer = questdlg(['There are several files with the naming pattern: ''' orig_pattern{1} '''. Should we merge them together ?'], 'Merging multiple files ?');
        do_merge = strncmp(answer,'Yes',3);
        if (~do_merge)
          merge_cmd = '';
          file_pattern = '';
        end
      end
    end
  end

  % Something went terribly wrong...
  if (isempty(format) | isempty(is_rgb))
    cd(curdir);
    error('Tracking:lociFormat', ['The metadata does not present the expected information: ''file format'' and ''RGB'' :\n\n' metadata]);
  end

  % Get the information out of the search results
  format = format{1}{1};
  is_rgb = strncmp(is_rgb{1}{1}, 'true', 4);
  is_grayscale = (~isempty(is_grayscale));

  if (strncmpi(format,'OME-TIFF',8) && isempty(file_pattern))

    % If it's already an OME-TIFF file, we do not need to convert it, so stop here
    newfile = fname;
    cd(curdir);
    newfile = relativepath(newfile);

    return;
  end

  % RGB need to be split up
  split_cmd = '';
  split_ext = '';
  if (is_rgb)
    split_cmd = '-separate ';
    split_ext = '_%c';

    if (is_grayscale)
      split_cmd = [split_cmd '-channel 0 '];
    end
  end
  nchan = (3-2*is_grayscale);

  % We create an OME-TIFF file
  if (use_tmp_folder)
    newname = fullfile(orig_path, [orig_name split_ext '.ome.tiff']);
  elseif (isempty(file_pattern))
    newname = [name split_ext '.ome.tiff'];
  else
    [file_path, file_name, file_ext] = fileparts(file_pattern{1});

    % If we merge the files, use also the folder name
    if (do_merge)
      [junk, tmp_name, junk] = fileparts(file_path);
      file_name = [tmp_name '_' file_name];
    end
    newname = fullfile(file_path, [file_name split_ext '.ome.tiff']);
  end

  % If the file already exists, we ask what to do
  if(exist(strrep(newname, '%c', '0'),'file'))

    % We initially do not know what to do
    answer = 0;
    if (forced)
      answer = 2;
    end

    % Creat the fancy name for display (otherwise it thinks they are LaTeX commands)
    [junk, tmp_name, junk] = fileparts(strrep(newname, split_ext, ''));
    printname = strrep(tmp_name(1:end-4),'_','\_');

    % We do not accept "empty" answers
    while (answer == 0)
      answer = menu(['The OME-TIFF version of ' printname ' already exists, overwrite it ?'],'Yes','No');
    end

    % Act accorindly to the policy
    switch answer

      % Delete the current files (did not dare overwriting it directly)
      case 1
        delete(strrep(newname, '%c', '*'));

      % Otherwise we can stop here
      case 2
        cd(curdir);

        % Store the new name(s)
        if (is_rgb)
          newfile = cell(nchan, 1);

          for c=0:nchan-1
            newfile{c+1} = relativepath(strrep(newname, '%c', num2str(c)));
          end
        else
          newfile = relativepath(newname);
        end

        return;
    end
  end

  % This also takes quite some time, so warn the user
  hInfo = warndlg('Converting to OME-TIFF, please wait...', 'Converting movie...');

  % Call directly the command line tool to do the job
  if (ispc)
    cmd_newname = ['"' newname '"'];
    [res, infos] = system(['bfconvert.bat ' merge_cmd split_cmd cmd_name ' ' cmd_newname]);
  else
    cmd_newname = strrep(newname,' ','\ ');
    [res, infos] = system(['./bfconvert ' merge_cmd split_cmd cmd_name ' ' cmd_newname]);
  end

  % Delete the temporary folder
  if (use_tmp_folder)
    delete(fullfile(tmp_folder, '*'));
    rmdir(tmp_folder);
  end

  % Delete the window if need be
  if (ishandle(hInfo))
    delete(hInfo);
  end

  % Check if an error occured
  if (res ~= 0)
    cd(curdir);
    error(infos);
  end

  % Store the new name in relative path and come back to the original folder
  cd(curdir);

  % Store the new name(s)
  if (is_rgb)
    newfile = cell(nchan, 1);

    for c=0:nchan
      newfile{c+1} = relativepath(strrep(newname, '%c', num2str(c)));
    end
  else
    newfile = relativepath(newname);
  end

  return;
end

function metadata = find_metadata(filename, metadata)
% This function tries to identify more suitable metadata. For now
% on, the following metadata are supported:
%   - Leica Application Suite ".las"
%   - Leica Application Suite ".cal.xml"
%   - uManager "metadata.txt"
%   - files manually placed in the "Metadata" folder and named
%     as the recording

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
  end

  return;
end
