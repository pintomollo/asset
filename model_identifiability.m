function [chi2ple, psple, errors] = model_identifiability(param_set, temp, nsteps)

  if (nargin == 0)
    param_set = 1;
    temp = 10;
    nsteps = 20;
  elseif (nargin == 1)
    temp = 10;
    nsteps = 20;
  elseif (nargin == 2)
    nsteps = 20;
  end

  test_distribution = false;
  if (param_set < 0)
    param_set = -param_set;
    test_distribution = true;
  end

  nestim = 200;

  switch param_set
    case 2
      fit_params = [4 5 10 11];
    case 3
      fit_params = [4 5 6 10 11];
    case 4
      fit_params = [2 4 5 6 8 10 11];
    case 5
      fit_params = [4 5 6 10 11 12];
    case 6
      fit_params = [1:12];
    otherwise
      fit_params = [4 10];
  end

  opts = get_struct('modeling');
  opts = load_parameters(opts, 'goehring.txt');

  x0 = opts.init_params;
  x0 = repmat(x0, [opts.nparticles, 1]);

  flow = opts.advection_params;
  if (size(flow, 1) ~= size(x0, 1))
    [X, Y] = meshgrid([1:size(flow, 1)], 1+([0:size(x0, 1)-1]*(size(flow, 1)-1)/(size(x0, 1)-1)).');
    flow = bilinear_mex(flow, X, Y, [2 2]);
  end

  ml_params = [opts.diffusion_params; ...
                opts.reaction_params(1:end-2, :)];

  [orig, orig_t] = simulate_model(x0, [ml_params; opts.reaction_params(end-1:end, :)], opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow, opts.user_data);
  orig = orig((end/2)+1:end, :);

  penalty = -((3*max(orig(:)))^2)*opts.nparticles/2;

  rescaling = 10.^(floor(log10(ml_params)));
  rescaling(ml_params == 0) = 1;
  ml_params = ml_params ./ rescaling;

  uuid = now + cputime;
  RandStream.setDefaultStream(RandStream('mt19937ar','Seed',uuid));
  
  if (test_distribution)
    opt = optimset('Display','off', 'Algorithm', 'levenberg-marquardt', 'MaxFunEvals', 5000, 'MaxIter', 5000);
    chi2diff = NaN(nestim, 1);
    noiseless = orig;
    err_count = 0;
    func_evals = 0;
    tmp_params = ml_params;

    for i=1:nestim
      orig = noiseless .* (1+0.2*randn(size(noiseless)));
      noisy_params = ml_params(fit_params) .* (1+0.1*randn(size(fit_params)));

      [best, chi2] = lsqnonlin(@chi2score, noisy_params, [], [], opt);
      L = sum(((orig(:) - noiseless(:)) / temp).^2);
      chi2diff(i) = L - chi2;
      fprintf(1, '.');
    end
    
    save(['ple-' num2str(uuid) '.mat'], 'chi2diff', 'func_evals', 'err_count', 'temp', 'fit_params');
  else

    nparams = numel(fit_params);
    chi2ple = cell(nparams, 1);
    psple = cell(nparams, 1);
    errors = NaN(nparams, 2);

    for i=1:nparams
      tmp_params = ml_params;
      err_count = 0;
      func_evals = 0;
      [chi2ple{i}, psple{i}] = ple(ml_params(fit_params), i, nsteps, @chi2score);
      errors(i, :) = [err_count func_evals];
    end

    save(['ple-' num2str(uuid) '.mat'], 'chi2ple', 'psple', 'errors', 'temp', 'fit_params');
  end

  return;

  function L = chi2score(params)

    func_evals = func_evals + 1;

    tmp_params(fit_params) = params;
    [res, t] = simulate_model(x0, [tmp_params .* rescaling; opts.reaction_params(end-1:end, :)], opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow, opts.user_data);
    res = res((end/2)+1:end, :);
    if (length(t) ~= length(orig_t))
      res = interp1q(t.', res.', orig_t.').';
    end

    L = sum(((orig - res) / temp).^2);
    err_count = err_count + any(isnan(L(:)));
    L(isnan(L)) = penalty;
    L = sum(L);

    return;
  end
end
