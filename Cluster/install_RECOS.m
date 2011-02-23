function install_RECOS(CPUUID)

  addpath('./CelegansAnalysis/');
  addpath('./CelegansAnalysis/libraries');
  addpath('./CelegansAnalysis/helpers');

  opts = get_struct('RECOS',1);
  opts.verbosity = 0;
  opts.auto_save = false;

  if (nargin == 1)
    opts.uuid = CPUUID;

    if (ischar(opts.uuid))
      opts.uuid = str2double(opts.uuid);
    end
  end

  RandStream.setDefaultStream(RandStream('mt19937ar','Seed',opts.uuid));

  if (exist('TmpData/tmpmat16.tmp', 'file'))
    opts.segmentation_type = 'dic';

    if (rand(1) > 0.5)
      opts.ml_type = 'cortex';
    else
      opts.ml_type = 'eggshell';
    end
  elseif (exist('TmpData/tmpmat17.tmp', 'file'))
    opts.segmentation_type = 'markers';
    opts.ml_type = 'cortex';
  elseif (exist('TmpData/tmpmat18.tmp', 'file'))
    opts.segmentation_type = 'markers';
    opts.ml_type = 'eggshell';
  else
    disp('No data file....');
    quit;
  end

  t = floor(rand(1) * 3);

  switch (t)
    case 0
      opts.do_ml = 'cmaes';
    case 1
      opts.do_ml = 'godlike';
    case 2
      opts.do_ml = 'pso';
  end

  disp(['ML: ' opts.do_ml '; Img: ' opts.segmentation_type '; Type: ' opts.ml_type]);
  
  RECOS('GZ920_MERGE_RECOS-', opts, 'config_file','optimized_params');
  quit

  return;
end
