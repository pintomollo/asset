function install_ASSET

  current_dir = '';
  if (numel(dir('get_struct.m')) == 1)
    current_dir = pwd;
    cd ..;
  else
    current_dir = 'celegans-analysis';
  end

  if (exist('get_struct') ~= 2)
    if (~isempty(current_dir))
      cd(current_dir);
    end

    here = pwd;
    addpath(here);
    addpath(fullfile(here, 'libraries'));
    addpath(fullfile(here, 'helpers'));

    cd ..;
  end

  if (exist('bilinear_mex') ~= 3)
    if (~isempty(current_dir))
      cd(current_dir);
    end
    
    try
      mex -setup;
      if (ispc)
        mex MEX\bilinear_mex.c;
        mex MEX\dp_score_mex.c;
        mex MEX\gaussian_mex.c MEX\gaussian_smooth.c;
        mex MEX\median_mex.c MEX\ctmf.c;
        mex MEX\imadm_mex.c MEX\gaussian_smooth.c;
      else
        mex MEX/bilinear_mex.c;
        mex MEX/dp_score_mex.c;
        mex MEX/gaussian_mex.c MEX/gaussian_smooth.c;
        mex MEX/median_mex.c MEX/ctmf.c;
        mex MEX/imadm_mex.c MEX/gaussian_smooth.c;
      end
    catch ME
      warning('Could not compile the MEX funtions, ASSET can still run without them but more slowly. Consider fixing this problem for more efficiency.');
    end

    cd ..;
  end

  try
    loci.formats.ImageReader();
    loci.formats.MetadataTools.createOMEXMLMetadata();
  catch ME
    warning('LOCI Bio-formats library does not seem to be installed, this will impair correct conversion of the images.\nVisit http://www.loci.wisc.edu/software/bio-formats and install the JAR libraries.')
  end

  if (~exist('TmpData', 'dir'))
    mkdir('TmpData');
  end

  return;
end
