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
  uuids = cell(fitting.nfits, 1);

  if (strncmp(fitting.type, 'data', 4))
    if (isempty(fitting.ground_truth))
      error('Data is missing for the fitting !');
    end

    embryo_size = range(fitting.x_pos);
    opts.reaction_params(end, :) = embryo_size;
    opts.boundaries = [0 embryo_size];
    opts.x_step = diff(opts.boundaries)/(opts.nparticles-1);
  end

  restart = false;
  if (~isempty(opts.init_func) & isa(opts.init_func, 'function_handle'))
    x0 = opts.init_func(opts);
    restart = true;
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
    nmaintenance = min(size(fitting.ground_truth, 2), 50) - 1;
    fitting.ground_truth = fitting.ground_truth(:, end-nmaintenance:end);
    fitting.t_pos = fitting.t_pos(end-nmaintenance:end);
    fitting.aligning_type = 'end';
    fitting.fit_flow = false;
  end

  size_data = size(fitting.ground_truth);

  if (strncmp(fitting.aligning_type, 'domain', 6))
    opts_expansion = load_parameters(opts, 'domain_expansion.txt');
    f = domain_expansion((fitting.ground_truth(1:end/2, :) + fitting.ground_truth((end/2)+1:end, :)).', size_data(1)/2, size_data(2), opts_expansion);
    frac_indx = find(f > fitting.fraction, 1, 'first');
  end

  linear_truth = fitting.ground_truth(:);
  simul_pos = ([0:opts.nparticles-1] * opts.x_step).';
  penalty = ((median(linear_truth))^2)*opts.nparticles;
  ndata = length(fitting.x_pos);

  full_error = penalty * size_data(2) * 10;
  log_error = -log(penalty);

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
    range_data = fitting.data_noise * range(linear_truth);
  end

  for f = 1:fitting.nfits
    uuids{f} = num2str(now + cputime);
    opt.LogFilenamePrefix = ['adr-kymo-' uuids{f} '_'];

    tmp_params = ml_params;

    p0 = ml_params(fit_params);

    if (fitting.fit_flow)
      p0 = [p0 1];
    else
      flow_scale = 1;
    end

    p0 = p0 .* (1+fitting.init_noise*randn(size(p0)));
    nparams = length(p0);

    if (strncmp(fitting.type, 'simulation', 5))
      fitting.ground_truth = noiseless + range_data*randn(size_data);
    end

    if (fitting.fit_flow)
      display(['Fitting ' num2str(nparams) ' parameters (' num2str(fit_params) ' & flow):']);
    else
      display(['Fitting ' num2str(nparams) ' parameters (' num2str(fit_params) '):']);
    end
    
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
      correct = true;

      if (fitting.fit_flow)
        flow_scale = p_all(end, i);
        tmp_params(fit_params) = p_all(1:end-1, i);
      else
        tmp_params(fit_params) = p_all(:, i);
      end

      if (restart)
        opts.reaction_params = tmp_params(2:end, :);
        [x0, correct] = opts.init_func(opts);
      end

      [res, t] = simulate_model(x0, tmp_params .* rescaling, opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow * flow_scale, opts.user_data, opts.max_iter);
      res = res((end/2)+1:end, :);
      if (opts.nparticles ~= ndata)
        res = interp1q(simul_pos, res, fitting.x_pos.');
      end

      res = [res; res];

      switch fitting.aligning_type
        case 'best'
          if (length(t) < size_data(2))
            cc = normxcorr2(res, fitting.ground_truth); 

            %[max_cc, imax] = max(cc(size_data(1), floor(size(res, 2)/2)+1:end));
            %corr_offset = [(imax-floor(size_data(2)/2))];
            
            [max_cc, imax] = max(cc(size(res, 1), :));
            corr_offset = -(imax-size(res, 2));
          elseif (all(isfinite(res)))
            cc = normxcorr2(fitting.ground_truth, res); 

            %[max_cc, imax] = max(cc(size_data(1), floor(size_data(2)/2)+1:end));
            %corr_offset = -[(imax-floor(size(res, 2)/2))];
            [max_cc, imax] = max(cc(size_data(1), :));
            corr_offset = -(imax-size_data(2));
          else
            corr_offset = 0;
          end
        case 'domain'
          f = domain_expansion(res(1:end/2, :).', size(res, 1)/2, size(res,2), opts_expansion);

          if (isnan(f(end)))
            err_all(i) = Inf;
            continue;
          end
          findx = find(f > fitting.fraction, 1, 'first');
          corr_offset = frac_indx - findx + 1;
          
          if (isempty(corr_offset))
            corr_offset = NaN;
          end
        case 'end'
          corr_offset = size_data(2) - size(res, 2) + 1;
      end

      gindxs = [0:size(res, 2)-1];
      goods = ((gindxs + corr_offset) > 0 & (gindxs + corr_offset) <= size_data(2));

      if (sum(goods) < 10)
        err_all(i) = Inf;
      else
        tmp = ones(size_data)*res(1, 2);
        tmp(:, corr_offset+gindxs(goods)) = res(:, gindxs(goods) + 1);

        res = tmp;

        if (fitting.scale_data)
          gres = ~isnan(res(:));
          c = [ones(sum(gres), 1), res(gres)] \ linear_truth(gres);
          
          if (c(2) <= 0)
            err_all(i) = Inf;
            continue;
          end

          res = c(1) + c(2)*res;
          
          %{
          %% Stop doing the full robust fit, too slow and not so useful...
          try
          c = robustfit(res(:), linear_truth);

          if (c(2) <= 0)
            err_all(i) = Inf;
            continue;
          end

          res = c(1) + c(2)*res;
          catch
            %Aaaa
          end
          %}
        end
        
        tmp_err = sum((fitting.ground_truth - res).^2);
        tmp_err(~isfinite(tmp_err)) = penalty;
        err_all(i) = sum(tmp_err);

        err_all(i) = err_all(i) + penalty*correct;
        
        %imagesc([res; fitting.ground_truth]);title([num2str(err_all(i)) ' (' num2str(log10(err_all(i))) ')']);drawnow

        bads = (tmp_params < 0);

        if (any(bads))
          err_all(i) = err_all(i) + sum(sum(exp(log_error*tmp_params(bads)), 2), 1);
        end
      end
    end

    if (flip)
      err_all = err_all.';
    end

    if (all(isinf(err_all)))
      % ML crashes when all values are non-numerical
      err_all(1) = full_error;
    end

    return;
  end
end
