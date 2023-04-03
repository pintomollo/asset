function install_ASSET(recompile)
  % INSTALL_ASSET adds the package to the work path, creates the required directories,
  % handles the dependencies and compiles the necessary MEX libraries.
  %
  % INSTALL_ASSET('UNINSTALL') removes this package from the work path.
  %
  % INSTALL_ASSET(RECOMPILE) if RECOMPILE is TRUE, recompiles all MEX files.
  %
  % Blanchoud group, University of Fribourg
  % Simon Blanchoud
  % 27.03.2023

  % By default, do not recompile everything
  if (nargin == 0)
    recompile = false;

    % In that particular case, delete the directories
  elseif (strcmp(recompile, 'uninstall'))
    cell_folder = which('install_ASSET');
    if (~isempty(cell_folder))
      [current_dir, junk, junk] = fileparts(cell_folder);
      folders = strsplit(path(), ':');

      % Remove the used directories from the work path
      for i=1:length(folders)
        if (strncmp(folders{i}, current_dir, length(current_dir)))
          rmpath(folders{i});
        end
      end
      savepath;

      % Remove the added java libraries if any
      fname = tilde_expand('~/.octaverc');
      fnew = tilde_expand('~/octave.rc');

      % We open the file and a new new
      fid = fopen(fname, 'rt');
      if fid >= 0
        fnw = fopen(fnew, 'wt');
        if fnw >= 0

          javastring = ['javaaddpath(''' current_dir];

          % And we copy everything that isn't in our directory
          line = fgetl(fid);
          while ischar(line)
            if (~strncmp(line, javastring, length(javastring)))
              fwrite(fnw, line);
            end
            line = fgetl(fid);
          end

          fclose(fnw);
          movefile(fnew, fname);
        end

        fclose(fid);
      end

      % Move to the root of SiBoRG
      [root_dir, junk, junk] = fileparts(current_dir);
      cd(root_dir);

      % Actually delete this package if required
      btn = questdlg ("Do you want to delete ASSET completely (cannot be undone)?", "DELETE FILES??", "Yes", "No", "No");
      if (strcmp (btn, "Yes"))
        disp('(not really) DELETED');
        % rmdir(current_dir, 's');
      end
    end

    return;
  end

  disp('Installing ASSET:');
  warning('off','all');

  % Start by moving inside the containing folder
  cell_folder = which('install_ASSET');
  [current_dir, junk, junk] = fileparts(cell_folder);
  [root_dir, junk, junk] = fileparts(current_dir);
  cd(current_dir);

  printf(' - updating octave path');
  % Add the proper directories to the work path
  subdirs = ls();
  folders = strsplit(path(), ':');
  for i=1:size(subdirs, 1)
    dir_name = strtrim(subdirs(i,:));
    if (isfolder(dir_name))
      if (dir_name(1) ~= '.' && dir_name(1) ~= '_')
        full_path = fullfile(current_dir, dir_name);
        if (~any(strncmp(folders, full_path, length(full_path))))
          addpath(full_path);
        end
      end
    end
    printf('.');
  end
  addpath(current_dir);
  savepath;
  printf('.done\n');

  printf(' - updating java path');

  % Added the required java libraries if any to octaverc
  fname = tilde_expand('~/.octaverc');

  % We open the file to add our files
  fid = fopen(fname, 'at');
  if fid >= 0

    % And we copy everything that is a jar in our libraries
    libraries = dir(fullfile('libraries', '*.jar'));
    fdisp(fid, ' ');

    % We both run the command and store it for the next boot
    for i=1:length(libraries)
      library = fullfile(pwd, 'libraries', libraries(i).name);
      fdisp(fid, ['javaaddpath(''' library ''');']);
      javaaddpath(library);
    end

    fclose(fid);
    printf('.');
  end
  printf('.done\n');

  printf(' - installing required packages');

  % Install the required packages
  if isempty(pkg('list', 'image'))
    pkg install -forge image;
  end
  printf('.');
  if isempty(pkg('list', 'io'))
    pkg install -forge io;
  end
  printf('.');
  if isempty(pkg('list', 'statistics'))
    pkg install -forge statistics;
  end
  printf('.');
  if isempty(pkg('list', 'control'))
    pkg install -forge control;
  end
  printf('.');
  if isempty(pkg('list', 'signal'))
    pkg install -forge signal;
  end
  printf('.');
  if isempty(pkg('list', 'tisean'))
    pkg install -forge tisean;
  end
  printf('.');
  if isempty(pkg('list', 'matgeom'))
    pkg install -forge matgeom;
  end

  %fid = fopen(fname, 'at');
  %if fid >= 0
  %  pkgs = pkg('list');

    %for i=1:length(pkgs)
    %  fdisp(fid, ['pkg load ' pkgs{i}.name ';']);
    %  pkg('load', pkgs{i}.name);
    %end

    %fclose(fid);
  %end
  printf('.done\n');

  printf(' - compiling MEX files');
  % Try to compile the necessary MEX files
  cd('MEX');
  files = ls('*_mex.c*');
  troubles = false;
  for i=1:size(files, 1)
    fname = strtrim(files(i,:));
    [junk, no_ext, ext] = fileparts(fname);
    if (recompile || exist(no_ext) ~= 3)
      failure = mex(fname);

      if (failure == 1)
        troubles = true;
        wtharning('SiBoRG:MEX', {'Could not compile the required MEX function!' ME.message});
      end
    end
    printf('.');
  end
  printf('.done\n');

  cd(root_dir);

  % These folders are required as well
  if (~exist('TmpData', 'dir'))
    mkdir('TmpData');
  end
  if (~exist('export', 'dir'))
    mkdir('export');
  end

  warning('on','all');
  warning('off', 'Octave:language-extension');
  warning('off', 'Octave:mixed-string-concat');
  warning('off', 'Octave:array-as-logical');
  warning('off', 'Octave:missing-semicolon');
  warning('off', 'Octave:shadowed-function');

  % Confirm to the user that everything went fine
  if (troubles)
    disp('Installation (almost) successful...');
  else
    disp('Installation successful !');
  end

  % Gnu GPL notice
  printf(['\nSimon''s Botrylloides Regeneration Group plateform,  Copyright (C) 2019  Simon Blanchoud\n', ...
    'This program comes with ABSOLUTELY NO WARRANTY;\n', ...
    'This is free software, and you are welcome to redistribute it\n', ...
    'under certain conditions; read licence.txt for details.\n\n']);

  return;
end
