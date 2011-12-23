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
    cd ..;
  else
    current_dir = 'celegans-analysis';
  end

  % If we can see this function, ASSET is already in the matlabpath
  if (exist('get_struct') ~= 2)
    if (~isempty(current_dir))
      cd(current_dir);
    end

    here = pwd;
    addpath(here);
    addpath(fullfile(here, 'libraries'));
    addpath(fullfile(here, 'helpers'));
    savepath

    cd ..;
  end

  % If this function exists, the MEX files were properly compiled already
  if (exist('bilinear_mex') ~= 3)
    if (~isempty(current_dir))
      cd(current_dir);
    end
    
    % Try to compile all the MEX functions
    try
      mex -setup;
      eval(['mex MEX' filesep 'bilinear_mex.c']);
      eval(['mex MEX' filesep 'dp_score_mex.c']);
      eval(['mex MEX' filesep 'gaussian_mex.c']);
      eval(['mex MEX' filesep 'median_mex.c']);
      eval(['mex MEX' filesep 'imadm_mex.c']);
      eval(['mex MEX' filesep 'emdc.c']);

      include_path = uigetdir('/usr/local/include', 'Select the directory containing the LEVMAR include files (.h)');

      eval(['mex -I' include_path ' -llapack -lblas -llevmar MEX' filesep 'fit_front_mex.c']);
    catch ME
      warning('Could not compile the MEX funtions, ASSET can still run without them but more slowly. Consider fixing this problem for more efficiency.');
    end

    cd ..;
  end

  % Try to use the java libraries to see if they are properly installed
  try
    loci.formats.ImageReader();
    loci.formats.MetadataTools.createOMEXMLMetadata();
  catch ME
    warning('LOCI Bio-formats library does not seem to be installed, this will impair correct conversion of the images.\nVisit http://www.loci.wisc.edu/software/bio-formats and install the JAR libraries.')
  end

  % This folder is required as well
  if (~exist('TmpData', 'dir'))
    mkdir('TmpData');
  end

  if (~exist('Config', 'dir'))
    cd('celegans-analysis');
    if (exist('Config', 'dir'))
      movefile('Config', '..');
    end
    cd ..;
  end

  return;
end
