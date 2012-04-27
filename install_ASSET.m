function install_ASSET
% INSTALL_ASSET adds the required directories to the matlabpath, compiles the MEX
% functions and checks for the existence of the required java libraries.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 19.05.2011

  % Start by checking where we are, we need to be inside the celegans-analysis folder
  current_dir = '';
  if (numel(dir('get_struct.m')) == 1)
    current_dir = pwd;
  else
    current_dir = 'celegans-analysis';
    cd(current_dir);
  end

  % If we can see this function, ASSET is already in the matlabpath
  if (exist('get_struct') ~= 2)
    here = pwd;
    addpath(here);
    addpath(fullfile(here, 'libraries'));
    addpath(fullfile(here, 'helpers'));
    savepath;
  end

  % If this function exists, the MEX files were properly compiled already
  is_first = true;
  had_problem = false;

  % Try to compile all the MEX functions
  if (exist('bilinear_mex') ~= 3)
    try
      if (is_first)
        mex -setup;
        is_first = false;
      end
      eval(['mex MEX' filesep 'bilinear_mex.c']);
    catch 
      had_problem = true;
    end
  end
  if (exist('dp_score_mex') ~= 3)
    try
      if (is_first)
        mex -setup;
        is_first = false;
      end
      eval(['mex MEX' filesep 'dp_score_mex.c']);
    catch 
      had_problem = true;
    end
  end
  if (exist('gaussian_mex') ~= 3)
    try
      if (is_first)
        mex -setup;
        is_first = false;
      end
      eval(['mex MEX' filesep 'gaussian_mex.c']);
    catch 
      had_problem = true;
    end
  end
  if (exist('median_mex') ~= 3)
    try
      if (is_first)
        mex -setup;
        is_first = false;
      end
      eval(['mex MEX' filesep 'median_mex.c']);
    catch 
      had_problem = true;
    end
  end
  if (exist('imadm_mex') ~= 3)
    try
      if (is_first)
        mex -setup;
        is_first = false;
      end
      eval(['mex MEX' filesep 'imadm_mex.c']);
    catch 
      had_problem = true;
    end
  end
  if (exist('emdc') ~= 3)
    try
      if (is_first)
        mex -setup;
        is_first = false;
      end
      eval(['mex MEX' filesep 'emdc.c']);
    catch 
      had_problem = true;
    end
  end
  if (exist('ellipse_distance_mex') ~= 3)
    try
      if (is_first)
        mex -setup;
        is_first = false;
      end
      eval(['mex MEX' filesep 'ellipse_distance_mex.c']);
    catch 
      had_problem = true;
    end
  end
  if (exist('symmetric_corrcoef_mex') ~= 3)
    try
      if (is_first)
        mex -setup;
        is_first = false;
      end
      eval(['mex MEX' filesep 'symmetric_corrcoef_mex.c']);
    catch 
      had_problem = true;
    end
  end
  if (exist('ellipse_ellipse_area_mex') ~= 3)
    try
      if (is_first)
        mex -setup;
        is_first = false;
      end
      eval(['mex MEX' filesep 'ellipse_ellipse_area_mex.c MEX' filesep 'ellipse_ellipse_overlap.c']);
    catch 
      had_problem = true;
    end
  end
  if (exist('fit_front_mex') ~= 3)
    try
      if (is_first)
        mex -setup;
        is_first = false;
      end
      include_path = uigetdir('/usr/local/include', 'Select the directory containing the LEVMAR include files (.h)');
      eval(['mex -I' include_path ' -llapack -lblas -llevmar MEX' filesep 'fit_front_mex.c']);
    catch 
      had_problem = true;
    end
  end
  
  if (had_problem)
    warning('Could not compile all the MEX funtions, ASSET can still run without them but more slowly. Consider fixing this problem for more efficiency.');
  end

  % Try to use the java libraries to see if they are properly installed
  if (exist('bfconvert.bat', 'file') ~= 2)
    if (exist('bftools', 'dir'))
      here = pwd;
      addpath(fullfile(here, 'bftools'));
      savepath;

      if (~exist('bfconvert.bat', 'file'))
        warning('ASSET:lociMissing', 'The LOCI command line tools are not present, this might be a problem if your recordings are not in TIFF format !\nYou can download them from http://www.loci.wisc.edu/bio-formats/downloads\nThen place the entire folder in the celegans-analysis folder.');
      end
    else
      warning('ASSET:lociMissing', 'The LOCI command line tools are not present, this might be a problem if your recordings are not in TIFF format !\nYou can download them from http://www.loci.wisc.edu/bio-formats/downloads\nThen place the entire folder in the celegans-analysis folder.');
    end
  end
  if (exist('bfconvert.bat', 'file') ~= 2)
    button = questdlg('Should ASSET try to install the LOCI command line tools ?');

    if (strncmpi(button, 'yes', 3))
      try
        rmdir('bftools', 's');
      catch
        % Nothing...
      end

      try
        unzip('http://loci.wisc.edu/files/software/bftools.zip', 'bftools');
        cd('bftools');
        urlwrite('http://loci.wisc.edu/files/software/loci_tools.jar', 'loci_tools.jar');
        cd ..;
        here = pwd;
        addpath(fullfile(here, 'bftools'));
        savepath;
      catch
        errs = lasterror;
        warning('ASSET:installLOCI', ['Installation failed for the following reason:\n' errs.message]);
      end

      if (exist('bfconvert.bat', 'file') == 2)
        msgbox('Installation successfull !');
      end
    end
  end

  % This folder is required as well
  if (~exist('TmpData', 'dir'))
    mkdir('TmpData');
  end

  cd ..;

  return;
end
