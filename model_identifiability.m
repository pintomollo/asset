function [chi2ple, psple, errors] = model_identifiability(param_set, temp, nsteps, noisy)

  if (nargin == 0)
    param_set = 1;
    temp = 10;
    nsteps = 20;
    noisy = false;
  elseif (nargin == 1)
    temp = 10;
    nsteps = 20;
    noisy = false;
  elseif (nargin == 2)
    nsteps = 20;
    noisy = false;
  elseif (nargin == 3)
    noisy = false;
  end

  test_distribution = false;
  if (param_set < 0)
    param_set = -param_set;
    test_distribution = true;
    noisy = true;
  end

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

  [orig, orig_t] = simulate_model(x0, [ml_params; opts.reaction_params(end-1:end, :)], opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow, opts.user_data, opts.max_iter);
  orig = orig((end/2)+1:end, :);

  if (noisy)
    size_data = size(orig);
    range_data = 0.15 * range(orig(:));
    %temp = range_data;
  end

  penalty = ((3*max(orig(:)))^2)*opts.nparticles/2;

  rescaling = 10.^(floor(log10(ml_params)));
  rescaling(ml_params == 0) = 1;
  ml_params = ml_params ./ rescaling;

  temp_norm = 1/(temp^2);

  uuid = now + cputime;
  rng(uuid, 'twister');
  
  if (test_distribution)
    noisy = false;

    opt = optimset('Display','off', 'Algorithm', 'levenberg-marquardt', 'MaxFunEvals', 5000, 'MaxIter', 5000);
    chi2diff = NaN(nsteps, 1);
    noiseless = orig;
    func_evals = 0;
    tmp_params = ml_params;
    err_count = 0;
    errors = NaN(nsteps, 1);
    
    for i=1:nsteps
      err_count = 0;
      orig = noiseless + range_data*randn(size_data);
      chi2 = chi2score(ml_params(fit_params));

      noisy_params = ml_params(fit_params) .* (1+0.1*randn(size(fit_params)));

      [best, junk, Loptim] = lsqnonlin(@chi2score, ml_params(fit_params), [], [], opt);
      %[best, junk, Loptim] = lsqnonlin(@chi2score, noisy_params, [], [], opt);
      %L = sum(((orig(:) - noiseless(:)) / temp).^2);
      %L = chi2score(noisy_params);
      chi2diff(i) = chi2 - Loptim;
      errors(i) = err_count;
      fprintf(1, '.');
    end
    
    save(['chi2-' num2str(uuid) '.mat'], 'chi2diff', 'func_evals', 'errors', 'temp', 'fit_params');
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
    [res, t] = simulate_model(x0, [tmp_params .* rescaling; opts.reaction_params(end-1:end, :)], opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow, opts.user_data, opts.max_iter);
    res = res((end/2)+1:end, :);
    if (length(t) ~= length(orig_t))
      res = interp1q(t.', res.', orig_t.').';
    end

    if (noisy)
      L = sum((orig + randn(size_data)*range_data - res).^2);
    else
      L = sum((orig - res).^2);
    end
    err_count = err_count + any(isnan(L(:)));
    L(~isfinite(L)) = penalty;
    L = sum(L) * temp_norm;

    return;
  end
end
