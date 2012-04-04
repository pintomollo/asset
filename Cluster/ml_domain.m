function ml_domain(CPUUID)

  addpath('./celegans-analysis/');
  addpath('./celegans-analysis/libraries');

  opts = get_struct('ASSET',1);
  opts.verbosity = 0;

  if (nargin == 1)
    opts.uuid = CPUUID;

    if (ischar(opts.uuid))
      opts.uuid = str2double(opts.uuid);
    end
  end

  RandStream.setDefaultStream(RandStream('mt19937ar','Seed',opts.uuid));

  %t = floor(rand(1) * 3);
  t = 1;

  switch (t)
    case 0
      opts.do_ml = 'cmaes';
    case 1
      opts.do_ml = 'godlike';
    case 2
      opts.do_ml = 'pso';
  end

  disp(['ML: ' opts.do_ml]);
  
  find_domain(opts);
  quit

  return;
end
