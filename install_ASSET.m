function install_ASSET(recompile)
% INSTALL_ASSET adds the required directories to the matlabpath, compiles the MEX
% functions and checks for the existence of the required java libraries.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 19.05.2011

  % By default, do not recompile everything
  if (nargin == 0)
    recompile = false;

  % In that particular case, delete the directories
  elseif (strcmp(recompile, 'uninstall'))
    cell_folder = which('install_ASSET');
    if (~isempty(cell_folder))
      [current_dir, junk, junk] = fileparts(cell_folder);
      [root_dir, junk, junk] = fileparts(current_dir);

      % Remove the used directories to MATLAB path
      rmpath(current_dir);
      rmpath(fullfile(current_dir, 'GUI'));
      rmpath(fullfile(current_dir, 'MEX'));
      rmpath(fullfile(current_dir, 'file_io'));
      rmpath(fullfile(current_dir, 'helpers'));
      rmpath(fullfile(current_dir, 'image_analysis'));
      rmpath(fullfile(current_dir, 'libraries'));
      rmpath(fullfile(current_dir, 'machine_learning'));
      rmpath(fullfile(current_dir, 'modules'));
      rmpath(fullfile(current_dir, 'pipeline'));
      rmpath(fullfile(current_dir, 'to_sort'));
      rmpath(fullfile(current_dir, 'unused'));
      rmpath(fullfile(current_dir, 'bftools'));
      savepath;

      return;
    end
  end

  % Start by moving inside the cell_tracking folder
  cell_folder = which('install_ASSET');
  [current_dir, junk, junk] = fileparts(cell_folder);
  [root_dir, junk, junk] = fileparts(current_dir);
  cd(current_dir);

  % Add the proper directories to MATLAB path
  addpath(current_dir);
  addpath(fullfile(current_dir, 'GUI'));
  addpath(fullfile(current_dir, 'MEX'));
  addpath(fullfile(current_dir, 'file_io'));
  addpath(fullfile(current_dir, 'helpers'));
  addpath(fullfile(current_dir, 'image_analysis'));
  addpath(fullfile(current_dir, 'libraries'));
  addpath(fullfile(current_dir, 'machine_learning'));
  addpath(fullfile(current_dir, 'modules'));
  addpath(fullfile(current_dir, 'pipeline'));
  addpath(fullfile(current_dir, 'to_sort'));
  savepath;

  % And for the LOCI as well
  if (exist(fullfile(current_dir, 'bftools'), 'dir'))
    addpath(fullfile(current_dir, 'bftools'));
    savepath;
  end

  % Otherwise, try to insall it !
  if (exist('bfconvert.bat', 'file') ~= 2)
    button = questdlg('Should we try to install the Bio-Formats command line tools ?');

    % Ask for the user to confirm this foolness
    if (strncmpi(button, 'yes', 3))
      try
        rmdir('bftools', 's');
      catch
        % Nothing...
      end

      % This looks like a permanent link... up to now at least
      try
        disp('Downloading Bio-Formats tools...')

        unzip('http://downloads.openmicroscopy.org/latest/bio-formats/artifacts/bftools.zip');
        addpath(fullfile(current_dir, 'bftools'));
        savepath;
      catch
        errs = lasterror;
        warning('ASSET:installLOCI', ['Installation failed for the following reason:\n' errs.message]);
      end

      % Amazing enough !!
      if (exist('bfconvert.bat', 'file') == 2)
        disp('Done !');
        disp(' ');
      else
        disp('Failed... try to get the Bio-Formats command-line tools from http://www.openmicroscopy.org and place it in the "cast" folder');
        disp(' ');
      end
    end
  end

  % Check if java is doing fine
  [res, info] = system('java -version');
  if isempty(strfind(info, 'java version "'))
    warning('ASSET:javaMissing', 'Java is either missing or not properly configured ! Make sure that its executable directory is included in your command path.');
  end

  % Check if the sparse 64bits flag is needed
  if ~isempty(strfind(computer(),'64'))
    mexopts = ' -largeArrayDims';
  else
    mexopts = '';
  end

  % Compile all files in MEX folder
  cd('MEX');
  [files] = get_mex_files('.', recompile);
  if (~isempty(files))
    problems = {};
    mex -setup;
    for i=1:length(files)
      try
        eval(['mex' mexopts ' ' files{i}]);
      catch ME
        problems{end+1} = ME.message;
      end
    end
    if (~isempty(problems))
      error('ASSET:install_ASSET', ['Could not compile the required MEX function(s)!\n' strjoin(problems, '\n')]);
    end
  end

  % Back to the original folder
  cd(root_dir);

  % This folder is required as well
  if (~exist('TmpData', 'dir'))
    mkdir('TmpData');
  end
  if (~exist('export', 'dir'))
    mkdir('export');
  end

  % Confirm to the user that everything went fine
  disp('Installation successful !');

  % Gnu GPL notice
  fprintf(1, ['\nASSET: ...,\n', ...
    'Copyright (C) 2014  Simon Blanchoud\n', ...
    'This program comes with ABSOLUTELY NO WARRANTY;\n', ...
    'This is free software, and you are welcome to redistribute it\n', ...
    'under certain conditions; read licence.txt for details.\n\n']);

  % First step !
  disp('Start using ASSET by calling "ASSET_GUI"');

  return;
end

% Get all the files to be compiled
function [files] = get_mex_files(dir_path, recompile)

  % All potential files
  candidates = [dir(fullfile(dir_path, '*.c')); dir(fullfile(dir_path, '*.cc')); dir(fullfile(dir_path, '*.cpp'))];
  files = {};
  ext = mexext();

  % Loop over them all
  for i=1:length(candidates)
    [junk, name, junk] = fileparts(candidates(i).name);

    % If there is no executable, then we might want to compile it
    if (~exist(fullfile(dir_path, [name '.' ext]), 'file') || recompile)
      curr_file = fullfile(dir_path, candidates(i).name);

      % Check if the file contains the mex main function
      content = fileread(curr_file);
      indx = strfind(content, 'void mexFunction(');

      % Store the valid files
      if (~isempty(indx))
        files{end+1} = curr_file;
      end
    end
  end

  return;
end
