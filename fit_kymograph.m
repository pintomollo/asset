function uuids = fit_kymograph(fitting, opts)

  if (nargin == 0)
    fitting = get_struct('fitting');
    opts = get_struct('modeling');
    opts = load_parameters(opts, 'goehring.txt');
    if (~fitting.fit_full)
      opts = load_parameters(opts, 'maintenance.txt');
    end
  elseif (nargin == 1)
    if (isfield(fitting, 'nparticles'))
      opts = fitting;
      fitting = get_struct('fitting');
    else
      opts = get_struct('modeling');
      opts = load_parameters(opts, 'goehring.txt');
      if (~fitting.fit_full)
        opts = load_parameters(opts, 'maintenance.txt');
      end
    end
  end

  switch fitting.parameter_set
    case 2
      fit_params = [4 5 12 13];
    case 3
      fit_params = [4 5 6 12 13];
    case 4
      fit_params = [2 4 5 6 10 12 13];
    case 5
      fit_params = [4 5 6 12 13 14];
    case 6
      fit_params = [1:16];
    otherwise
      fit_params = [4 12];
  end

  RandStream.setDefaultStream(RandStream('mt19937ar','Seed', now + cputime));
  uuids = NaN(fitting.nfits, 1);

  if (strncmp(fitting.type, 'data', 4))
    if (isempty(fitting.ground_truth))
      error('Data is missing for the fitting !');
    end

    embryo_size = 2*range(fitting.x_pos);
    opts.reaction_params(end, :) = embryo_size;
    opts.boundaries = [0 embryo_size/2];
    opts.x_step = diff(opts.boundaries)/(opts.nparticles-1);
  end

  if (~isempty(opts.init_func) & isa(opts.init_func, 'function_handle'))
    x0 = opts.init_func(opts);
  else
    x0 = opts.init_params;
    x0 = repmat(x0, [opts.nparticles, 1]);
  end

  flow = opts.advection_params;
  if (size(flow, 1) ~= size(x0, 1))
    [X, Y] = meshgrid([1:size(flow, 1)], 1+([0:size(x0, 1)-1]*(size(flow, 1)-1)/(size(x0, 1)-1)).');
    flow = bilinear_mex(flow, X, Y, [2 2]);
  end

  ml_params = [opts.diffusion_params; ...
                opts.reaction_params];

  if (strncmp(fitting.type, 'simulation', 5))
    [fitting.ground_truth, fitting.t_pos] = simulate_model(x0, ml_params, opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow, opts.user_data, opts.max_iter);
    fitting.ground_truth = fitting.ground_truth((end/2)+1:end, :);
    fitting.ground_truth = [fitting.ground_truth; fitting.ground_truth];
    fitting.x_pos = [0:opts.nparticles-1] * opts.x_step;
  end

  if (~fitting.fit_full)
    nmaintenance = min(size(fitting.ground_truth, 2), 100) - 1;
    fitting.ground_truth = fitting.ground_truth(:, end-nmaintenance:end);
    fitting.t_pos = fitting.t_pos(end-nmaintenance:end);
    fitting.aligning_type = 'end';
  end

  linear_truth = fitting.ground_truth(:);
  simul_pos = ([0:opts.nparticles-1] * opts.x_step).';
  penalty = ((3*max(linear_truth))^2)*opts.nparticles/2;
  ndata = length(fitting.x_pos);

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

  if (strncmp(fitting.type, 'simulation', 5))
    noiseless = fitting.ground_truth;
    size_data = size(noiseless);
    range_data = fitting.data_noise * range(linear_truth);
  end

  for f = 1:fitting.nfits
    uuids(f) = now + cputime;
    opt.LogFilenamePrefix = ['adr-kymo-' num2str(uuids(f)) '_'];

    tmp_params = ml_params;

    p0 = ml_params(fit_params);
    p0 = p0 .* (1+fitting.init_noise*randn(size(p0)));
    nparams = length(p0);

    if (strncmp(fitting.type, 'simulation', 5))
      fitting.ground_truth = noiseless + range_data*randn(size_data);
    end

    display(['Fitting ' num2str(nparams) ' parameters (' num2str(fit_params) '):']);
    
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
      [res, t] = simulate_model(x0, tmp_params .* rescaling, opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow, opts.user_data, opts.max_iter);
      res = res((end/2)+1:end, :);
      if (length(t) ~= length(fitting.t_pos))
        res = interp1q(t.', res.', fitting.t_pos.').';
      end
      if (opts.nparticles ~= ndata)
        res = interp1q(simul_pos, res, fitting.x_pos.');
      end

      %if (~fitting.fit_full)
      %  res = res(:, end-100:end);
      %end
      res = [res; res];

      switch fitting.aligning_type
        case 'best'
          %cc = normxcorr2(ground_truth, s); 
          %[max_cc, imax] = max(cc(:,size(ground_truth, 2)));
          %corr_offset = [(imax-size(ground_truth,1))];
        case 'end'
          %corr_offset = size(s, 1) - size(ground_truth, 1) + 1;
      end

      if (fitting.scale_data)
        try
        c = robustfit(res(:), linear_truth);

        res = c(1) + c(2)*res;
        catch
          %Aaaa
        end
      end

      tmp_err = sum((fitting.ground_truth - res).^2);
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
