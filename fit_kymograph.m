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
  else
    orig_opts = get_struct('modeling');
    orig_opts = load_parameters(orig_opts, 'goehring.txt');
    if (fitting.fit_relative)
      orig_opts.reaction_params(3,:) = orig_opts.reaction_params(3,:) .* (orig_opts.reaction_params(5,[2 1]).^(orig_opts.reaction_params(4,:))) ./ orig_opts.reaction_params(4,[2 1]);
      orig_opts.reaction_params(5,:) = 1;
    end
    orig_params = [orig_opts.diffusion_params; ...
                orig_opts.reaction_params];

    orig_scaling = 10.^(floor(log10(orig_params)));
    orig_scaling(orig_params == 0) = 1;
  end

  if (fitting.fit_relative)
    opts.reaction_params(3,:) = opts.reaction_params(3,:) .* (opts.reaction_params(5,[2 1]).^(opts.reaction_params(4,:))) ./ opts.reaction_params(4,[2 1]);
    opts.reaction_params(5,:) = 1;
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

  rng(now + cputime, 'twister');
  uuids = cell(fitting.nfits, 1);
  dp = 1e-3;

  %tail_coeff = 0.25;
  %stable_coeff = 2;
  %close_coeff = 10;

  normalization_done = false;

  if (strncmp(fitting.fitting_type, 'dram', 4))
    drscale  = 2; 
    adaptint = 100;
  end

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
    fitting.ground_truth = fitting.ground_truth(:, end-nmaintenance:end, :);
    fitting.t_pos = fitting.t_pos(end-nmaintenance:end);
    fitting.aligning_type = 'end';
    fitting.fit_flow = false;
  end

  size_data = size(fitting.ground_truth);

  if (strncmp(fitting.aligning_type, 'domain', 6))
    opts_expansion = load_parameters(get_struct('ASSET'), 'domain_expansion.txt');
    [f, frac_width, full_width] = domain_expansion(mymean(fitting.ground_truth(1:end/2, :, :), 3).', mymean(fitting.ground_truth((end/2)+1:end, :, :), 3).', size_data(1)/2, size_data(2), opts_expansion);
    frac_indx = find(f > fitting.fraction, 1, 'first');

    if (strncmp(fitting.scale_type, 'normalize', 10))
      fitting.ground_truth = normalize_domain(fitting.ground_truth, f*frac_width, opts_expansion);
      normalization_done = true;
    end
  end

  if (~normalization_done & strncmp(fitting.scale_type, 'normalize', 10))
    opts_expansion = load_parameters(get_struct('ASSET'), 'domain_expansion.txt');
    [f, frac_width, full_width] = domain_expansion(mymean(fitting.ground_truth(1:end/2, :, :), 3).', mymean(fitting.ground_truth((end/2)+1:end, :, :), 3).', size_data(1)/2, size_data(2), opts_expansion);
    if (~fitting.fit_full)
      f(:) = 1;
    end
    fitting.ground_truth = normalize_domain(fitting.ground_truth, f*frac_width, opts_expansion);
    normalization_done = true;
  end

  if (strncmp(fitting.aligning_type, 'lsr', 3))
    sindx = round(size(fitting.ground_truth, 2) / 2);
    mean_ground_truth = mymean(fitting.ground_truth, 3);
  end

  multi_data = (numel(size_data) > 2 & size_data(3) > 1);
  if (multi_data)
    nlayers = size_data(3);
    size_data = size_data(1:2);
  end

  linear_truth = fitting.ground_truth(:);
  linear_goods = isfinite(linear_truth);
  simul_pos = ([0:opts.nparticles-1] * opts.x_step).';
  penalty = ((max(linear_truth(linear_goods)))^2)*opts.nparticles;
  ndata = length(fitting.x_pos);

  full_error = penalty * size_data(2) * 10;
  log_error = -log(penalty);

  rescaling = orig_scaling;
  %rescaling = 10.^(floor(log10(ml_params)));
  %rescaling(ml_params == 0) = 1;
  ml_params = ml_params ./ rescaling;

  if (strncmp(fitting.type, 'simulation', 5))
    noiseless = fitting.ground_truth;
    range_data = fitting.data_noise * range(linear_truth);
  end

  for f = 1:fitting.nfits
    uuids{f} = num2str(now + cputime);
    log_name = ['adr-kymo-' uuids{f} '_'];

    tmp_params = ml_params;

    p0 = ml_params(fit_params);

    if (fitting.fit_flow)
      p0 = [p0 1];
    else
      flow_scale = 1;
    end

    p0 = p0 .* (1+fitting.init_noise*randn(size(p0))/sqrt(length(p0)));
    nparams = length(p0);

    if (strncmp(fitting.type, 'simulation', 5))
      fitting.ground_truth = noiseless + range_data*randn(size_data);
    end

    if (fitting.fit_flow)
      display(['Fitting ' num2str(nparams) ' parameters (' num2str(fit_params) ' & flow):']);
    else
      display(['Fitting ' num2str(nparams) ' parameters (' num2str(fit_params) '):']);
    end
    p0 = sqrt(p0(:));

    switch (fitting.fitting_type)
      case 'cmaes'

        opt = cmaes('defaults');
        opt.MaxFunEvals = fitting.max_iter;
        opt.TolFun = dp;
        opt.SaveFilename = '';
        opt.SaveVariables = 'off';
        opt.EvalParallel = 'yes';
        opt.LogPlot = 0;
        opt.LogFilenamePrefix = log_name;

        [p, fval, ncoutns, stopflag, out] = cmaes(@error_function, p0(:), 0.25, opt); 
      case 'pso'
        opt = [1 2000 24 0.5 0.5 0.7 0.2 1500 dp 250 NaN 0 0];
        opt(2) = fitting.max_iter;
        opt(9) = 1e-10;
        opt(12) = 1;
        opt(13) = 1;
     

        bounds = [zeros(nparams, 1) ones(nparams, 1)*10];
        [pbest, tr, te] = pso_Trelea_vectorized(@error_function, nparams, NaN, bounds, 0, opt, '', p0.', [log_name 'evol']);
        p = pbest(1:end-1);
      case 'godlike'
        opt = set_options('Display', 'on', 'MaxIters', fitting.max_iter, 'TolFun', dp, 'LogFile', [log_name 'evol']);
        which_algos = {'DE';'GA';'PSO';'ASA'};
        which_algos = which_algos(randperm(4));
        [p, fval] = GODLIKE(@error_function, 20, zeros(nparams,1), ones(nparams,1) * 10, which_algos, opt);

      case 'dram'
%        [estim_error, estim_sigma2] = error_function(p0(:));

        model.ssfun    = @error_function;

        params.par0    = p0(:); % initial parameter values
        params.n       = sum(linear_goods);
        params.n0      = params.n;  % prior for error variance sigma^2
%        params.sigma2  = estim_sigma2;  % prior accuracy for sigma^2

        options.nsimu    = fitting.max_iter;               % size of the chain
        options.adaptint = adaptint;            % adaptation interval
        options.drscale  = drscale;
        options.qcov     = fitting.step_size*eye(nparams).*2.4^2./nparams;      % initial proposal covariance 
        options.ndelays  = fitting.ndelays;
        options.log_file = [log_name 'evol'];

        % run the chain
        results = dramrun(model,[],params,options);
        p = results.mean;
      otherwise
        error [opts.do_ml ' machine learning algorithm is not implemented'];

        return;
    end

    p = p.^2;

    display(['Best (' num2str(p(:).') ')']);
  end

  return;

  function [err_all, sigma2] = error_function(varargin)

    sigma2 = 0;
    p_all = varargin{1};

    [curr_nparams, nevals] = size(p_all);
    flip = false;
    if (nevals == nparams)
      flip = true;
      p_all = p_all.';
      nevals = curr_nparams; 
    end
    err_all = NaN(1, nevals);

    p_all = p_all.^2;

    for i = 1:nevals
      correct = true;

      if (fitting.fit_flow)
        flow_scale = p_all(end, i);
        tmp_params(fit_params) = p_all(1:end-1, i);
      else
        tmp_params(fit_params) = p_all(:, i);
      end

      if (restart)
        opts.reaction_params = tmp_params(2:end, :) .* rescaling(2:end, :);
        [x0, correct] = opts.init_func(opts);

        if (~correct)
          err_all(i) = Inf;
          continue;
        end
      end

      normalization_done = false;
      if (fitting.fit_relative)
        [res, t] = simulate_model_rel(x0, tmp_params .* rescaling, opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow * flow_scale, opts.user_data, opts.max_iter);
      else
        [res, t] = simulate_model_mix(x0, tmp_params .* rescaling, opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow * flow_scale, opts.user_data, opts.max_iter);
      end
      
      %if (t(end) < opts.tmax)
      %  err_all(i) = Inf;
      %  continue;
      %end
      
      res = res((end/2)+1:end, :);
      if (opts.nparticles ~= ndata)
        res = interp1q(simul_pos, res, fitting.x_pos.');
      end
      
      res = [res; res];
      switch fitting.aligning_type
        case 'best'
          if (length(t) < size_data(2))
            cc = normxcorr2(res, mymean(fitting.ground_truth, 3)); 

            [max_cc, imax] = max(cc(size(res, 1), :));
            corr_offset = -(imax-size(res, 2));
          elseif (all(isfinite(res)))
            cc = normxcorr2(mymean(fitting.ground_truth, 3), res); 

            [max_cc, imax] = max(cc(size_data(1), :));
            corr_offset = -(imax-size_data(2));
          else
            corr_offset = 0;
          end
        case 'domain'
          if (size(res,2) <= 10)
            err_all(i) = Inf;
            continue;
          else
            [f, fwidth] = domain_expansion(res(1:end/2, :).', size(res, 1)/2, size(res,2), opts_expansion);

            if (strncmp(fitting.scale_type, 'normalize', 10))
              res = normalize_domain(res, f*fwidth, opts_expansion);
              normalization_done = true;
            end
          end

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
        case 'lsr'
          rindx = min(sindx, size(res, 2));

          [junk, junk2, rindx] = find_min_residue(mean_ground_truth, sindx, res, rindx, 0.95);
          corr_offset = sindx - rindx;
      end

      gindxs = [0:size(res, 2)-1];
      goods = ((gindxs + corr_offset) > 0 & (gindxs + corr_offset) <= size_data(2));

      if (sum(goods) < 10)
        err_all(i) = Inf;
      else
        if (~normalization_done & strncmp(fitting.scale_type, 'normalize', 10))
          [f, fwidth] = domain_expansion(res(1:end/2, :).', size(res, 1)/2, size(res,2), opts_expansion);
          if (~fitting.fit_full)
            f(:) = 1;
          end
          res = normalize_domain(res, f*fwidth, opts_expansion);
          normalization_done = true;
        end

        tmp = ones(size_data)*res(1, 2);
        tmp(:, corr_offset+gindxs(goods)) = res(:, gindxs(goods) + 1);

        res = tmp;
        if (multi_data)
          res = repmat(res, [1 1 nlayers]);
        end

        if (~normalization_done & strncmp(fitting.scale_type, 'best', 4))
          gres = ~isnan(res(:)) & linear_goods;
          c = [ones(sum(gres), 1), res(gres)] \ linear_truth(gres);
          
          if (c(2) <= 0)
            err_all(i) = Inf;
            continue;
          end

          res = c(1) + c(2)*res;
        end
        
        tmp_err = (fitting.ground_truth - res).^2;

        if (nargout == 2)
          tmp_mean = mymean(tmp_err(:));
          sigma2 = mymean((tmp_err(:) - tmp_mean).^2);
        end

        tmp_err(~linear_goods) = 0;
        err_all(i) = sum(tmp_err(:));

        %figure;
        %subplot(1,2,1);
        %imagesc(mymean(res, 3));
        %subplot(1,2,2);
        %imagesc(mymean(tmp_err, 3));
        %title([num2str(err_all(i)) ' ' num2str(p_all(:,i).')]);
        %keyboard
      end
    end

    err_all(isnan(err_all)) = Inf;

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

function domain = normalize_domain(domain, path, opts)

  path = path/opts.quantification.resolution;
  [h, w, f] = size(domain);
  h = h/2;
  pos_mat = repmat([1:h].', 1, w);
  mask = bsxfun(@le, pos_mat, path.');

  full_mask = repmat(flipud(mask), [2, 1, f]);
  bkg = mymean(domain(~full_mask));
  int = mymean(domain(full_mask));

  if (isempty(bkg))
    bkg = 0;
  end

  domain = (domain - bkg) / (int-bkg);

  return;
end
