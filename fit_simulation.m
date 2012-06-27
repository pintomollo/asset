function fit_simulation(param_set, init_noise, noise_data, nfits, opts)

  if (nargin < 5)
    opts = get_struct('modeling');
    opts = load_parameters(opts, 'goehring.txt');
  end

  scale_data = true;

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

  uuid = now + cputime;
  RandStream.setDefaultStream(RandStream('mt19937ar','Seed',uuid));

  x0 = opts.init_params;
  x0 = repmat(x0, [opts.nparticles, 1]);

  flow = opts.advection_params;
  if (size(flow, 1) ~= size(x0, 1))
    [X, Y] = meshgrid([1:size(flow, 1)], 1+([0:size(x0, 1)-1]*(size(flow, 1)-1)/(size(x0, 1)-1)).');
    flow = bilinear_mex(flow, X, Y, [2 2]);
  end

  ml_params = [opts.diffusion_params; ...
                opts.reaction_params(1:end-2, :)];

  [noiseless, orig_t] = simulate_model(x0, [ml_params; opts.reaction_params(end-1:end, :)], opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow, opts.user_data, opts.max_iter);
  %noiseless = noiseless((end/2)+1:end, :);
  noiseless = noiseless((end/2)+1:end, end-100:end);
  penalty = ((3*max(noiseless(:)))^2)*opts.nparticles/2;

  size_data = size(noiseless);
  range_data = noise_data * range(noiseless(:));

  rescaling = 10.^(floor(log10(ml_params)));
  rescaling(ml_params == 0) = 1;
  ml_params = ml_params ./ rescaling;

  opt = cmaes('defaults');
  opt.MaxFunEvals = 10000;
  opt.TolFun = 1e-5;
  opt.SaveFilename = '';
  opt.SaveVariables = 'off';
  opt.EvalParallel = 'yes';
  opt.LogPlot = 0;

  for f=1:nfits
    uuid = now + cputime;
    opt.LogFilenamePrefix = ['adr-fit-' num2str(uuid) '_'];

    tmp_params = ml_params;
    orig = noiseless + range_data*randn(size_data);

    p0 = ml_params(fit_params);
    p0 = p0 .* (1+init_noise*randn(size(p0)));
    nparams = length(p0);

    display(['Fitting ' num2str(nparams) ' parameters (' num2str(fit_params) '):']);
    display(['IC (' num2str(p0) ') : GT (' num2str(ml_params(fit_params)) ')']);
    
    [p, fval, ncoutns, stopflag, out] = cmaes(@error_function, p0(:), 0.5, opt);

    display(['Best (' num2str(p.') ')']);
  end

  return;

  function err_all = error_function(p_all)

    [curr_nparams, nevals] = size(p_all);
    flip = false;
    if (nevals == nparams)
      flip = true;
      p_all = p_all.';
      nevals = curr_nparams; 
    end
    err_all = NaN(1, nevals);

    for i = 1:nevals
      tmp_params(fit_params) = p_all(:, i);
      [res, t] = simulate_model(x0, [tmp_params .* rescaling; opts.reaction_params(end-1:end, :)], opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow, opts.user_data, opts.max_iter);
      res = res((end/2)+1:end, :);
      if (length(t) ~= length(orig_t))
        res = interp1q(t.', res.', orig_t.').';
      end
      res = res(:, end-100:end);

      if (scale_data)
        try
        c = robustfit(res(:), orig(:));

        res = c(1) + c(2)*res;
        catch
          %Aaaa
        end
      end

      tmp_err = sum((orig - res).^2);
      tmp_err(~isfinite(tmp_err)) = penalty;
      err_all(i) = sum(tmp_err);

      bads = (tmp_params < 0);

      if (any(bads))
        err_all(i) = err_all(i) + sum(sum(exp(-10*tmp_params(bads)), 2), 1);
      end
    end

    if (flip)
      err_all = err_all.';
    end

    return;
  end
end
