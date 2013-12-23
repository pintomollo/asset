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
    %orig_opts = get_struct('modeling');
    %orig_opts = load_parameters(orig_opts, 'goehring.txt');
    orig_opts = opts;
    if (fitting.fit_relative)
      orig_opts.reaction_params(3,:) = orig_opts.reaction_params(3,:) .* orig_opts.reaction_params(4,[2 1]);
      %orig_opts.reaction_params(3,:) = orig_opts.reaction_params(3,:) .* (orig_opts.reaction_params(5,[2 1]).^(orig_opts.reaction_params(4,:))) ./ orig_opts.reaction_params(4,[2 1]);
      %orig_opts.reaction_params(5,:) = 1;
    end
    orig_params = [orig_opts.diffusion_params; ...
                orig_opts.reaction_params];

    orig_scaling = 10.^(floor(log10(orig_params)));
    orig_scaling(orig_params == 0) = 1;
  end

  if (fitting.fit_relative)
    opts.reaction_params(3,:) = opts.reaction_params(3,:) ./ opts.reaction_params(4,[2 1]);
    %opts.reaction_params(3,:) = opts.reaction_params(3,:) .* (opts.reaction_params(5,[2 1]).^(opts.reaction_params(4,:))) ./ opts.reaction_params(4,[2 1]);
    %opts.reaction_params(5,:) = 1;
  end

  ngroups = length(fitting.ground_truth);

  all_params = false;

  noffsets = ngroups*strncmp(fitting.aligning_type, 'fitting', 7);
  [temperatures, junk, temp_indx] = unique(fitting.temperature);
  good_visc = (temperatures ~= opts.reaction_temperature);
  bad_flow = (temperatures == opts.flow_temperature);
  nvisc = sum(good_visc);
  ntemps = length(temperatures);
  viscosities = ones(1, ntemps);

  ml_params = [opts.diffusion_params; ...
                opts.reaction_params];

  [fit_params, fit_energy, fit_temperatures, fit_viscosity] = model_description(fitting.parameter_set);

  % In these case, we always need a model
  if (mod(fitting.parameter_set, 10) < 3)
    fitting.fit_model = true;
  end

  nrates = length(fit_params);
  nenergy = length(fit_energy);

  if (~fitting.fit_model && fit_temperatures)
    fit_energy = ones(1, nenergy*nvisc);
    fitting.fit_flow = any(~good_visc);
  end

  if (~isempty(fitting.init_pos))
    fitting.init_pos = fitting.init_pos(:).';

    if (iscell(fitting.init_pos))
      tmp_fit = fitting;
      uuids = [];

      for c=1:length(fitting.init_pos)
        tmp_fit.init_pos = fitting.init_pos{c};
        tmp_uuids = fit_kymograph(tmp_fit, opts);
        uuids = [uuids; tmp_uuids(:)];
      end

      return;
    elseif ((nrates + numel(fit_energy) + noffsets + fitting.fit_flow + fit_viscosity*nvisc) == numel(fitting.init_pos))

      fitting.init_pos(1:nrates) = fitting.init_pos(1:nrates) ./ orig_scaling(fit_params);
      fitting.init_pos(nrates+1:nrates+noffsets) = fitting.init_pos(nrates+1:nrates+noffsets) ./ fitting.offset_scaling;

      all_params = true;
    elseif ((nrates + numel(fit_energy) + noffsets + fit_viscosity*nvisc) == numel(fitting.init_pos))

      fitting.init_pos(1:nrates) = fitting.init_pos(1:nrates) ./ orig_scaling(fit_params);
      fitting.init_pos(nrates+1:nrates+noffsets) = fitting.init_pos(nrates+1:nrates+noffsets) ./ fitting.offset_scaling;
    elseif (nrates == numel(fitting.init_pos))
      fitting.init_pos = fitting.init_pos ./ orig_scaling(fit_params);
      fitting.init_pos = [fitting.init_pos fit_energy];
    elseif (4 + noffsets == numel(fitting.init_pos))
      warning('Assuming only the 4 main ones were provided !')
      tmp_params = ml_params;
      tmp_params([4 5 12 13]) = fitting.init_pos(1:4);
      tmp_params = tmp_params ./ orig_scaling;

      fitting.init_pos = [tmp_params(fit_params) fitting.init_pos(5:end)/fitting.offset_scaling fit_energy];
    elseif (4 == numel(fitting.init_pos))
      warning('Assuming only the 4 main ones were provided !')
      tmp_params = ml_params;
      tmp_params([4 5 12 13]) = fitting.init_pos(1:4);
      tmp_params = tmp_params ./ orig_scaling;

      fitting.init_pos = [tmp_params(fit_params) fit_energy];
    else
      warning('The provided initial position does not correspond to the dimensionality of the fit, ignoring it.');
      %keyboard
      fitting.init_pos = [];
    end
  end

  if (exist('rng'))
    rng(now + cputime, 'twister');
  else
    RandStream.setDefaultStream(RandStream('mt19937ar','seed',now + cputime));
  end
  uuids = cell(fitting.nfits, 1);

  %tail_coeff = 0.25;
  %stable_coeff = 2;
  %close_coeff = 10;

  normalization_done = false;

  if (strncmp(fitting.fitting_type, 'dram', 4))
    drscale  = 2; 
    %adaptint = 500;
    adaptint = 0;
  end

  if (~strncmp(fitting.type, 'simulation', 10))
    if (isempty(fitting.ground_truth))
      error('Data is missing for the fitting !');
    end

    %embryo_size = range(fitting.x_pos);
    %opts.reaction_params(end, :) = embryo_size;
    %opts.boundaries = [0 embryo_size];
    opts.boundaries = [0 opts.reaction_params(end,1)];
    opts.x_step = diff(opts.boundaries)/(opts.nparticles-1);
  end

  flow = opts.advection_params;
  if (size(flow, 1) ~= opts.nparticles)
    [X, Y] = meshgrid([1:size(flow, 2)], 1+([0:opts.nparticles-1]*(size(flow, 1)-1)/(opts.nparticles-1)).');
    flow = bilinear_mex(flow, X, Y, [2 2]);
  end

  if (fitting.display)
    tmp_h = findobj('Type', 'figure');
    hfig = tmp_h(1:min(end, 3));
  end

  for g=1:ngroups

    if (fitting.display)
      if (length(hfig) < g)
        hfig(g) = figure('Visible', 'off', 'Name', ['Group #' num2str(g)]);
      else
        set(hfig(g), 'Name', ['Group #' num2str(g)]);
      end
    end

    is_data = ~(strncmp(fitting.type, 'simulation', 10));
    if (~is_data)

        if (~isempty(opts.init_func) & isa(opts.init_func, 'function_handle'))
          x0 = opts.init_func(opts, fitting.fit_relative);
        else
          x0 = opts.init_params;
          x0 = repmat(x0, [opts.nparticles, 1]);
        end

        fitting.simulation_parameters = ml_params(fit_params);

        if (fitting.fit_relative)
           [fitting.ground_truth, fitting.t_pos] = simulate_model_real(x0, ml_params, opts.x_step, opts.tmax*0.75, opts.time_step, opts.output_rate, flow, opts.user_data, opts.max_iter);
        else
           [fitting.ground_truth, fitting.t_pos] = simulate_model_mix(x0, ml_params, opts.x_step, opts.tmax*0.75, opts.time_step, opts.output_rate, flow, opts.user_data, opts.max_iter);
        end

      %[fitting.ground_truth, fitting.t_pos] = simulate_model(x0, ml_params, opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow, opts.user_data, opts.max_iter);
      fitting.ground_truth = fitting.ground_truth((end/2)+1:end, :);
      fitting.ground_truth = {[fitting.ground_truth; fitting.ground_truth]};
      fitting.x_pos = {[0:opts.nparticles-1] * opts.x_step};

    elseif (~isempty(fitting.t_pos))
      if (numel(fitting.t_pos{g}) == 1)
        opts.output_rate(g) = fitting.t_pos{g};
      else
        opts.output_rate(g) = median(diff(fitting.t_pos{g}));
      end
    end

    if (~fitting.fit_full)
      nmaintenance = min(size(fitting.ground_truth{g}, 2), 50) - 1;
      fitting.ground_truth{g} = fitting.ground_truth{g}(:, end-nmaintenance:end, :);
      fitting.t_pos{g} = fitting.t_pos{g}(end-nmaintenance:end);
      fitting.aligning_type = 'end';
      fitting.fit_flow = false;
    end

    [n,m,o] = size(fitting.ground_truth{g});
    size_data(g,1:3) = [n,m,o];

    %if (strncmp(fitting.aligning_type, 'domain', 6))
      opts_expansion = load_parameters(get_struct('ASSET'), 'domain_expansion.txt');
      [f, frac_width, full_width] = domain_expansion(mymean(fitting.ground_truth{g}(1:end/2, :, :), 3).', mymean(fitting.ground_truth{g}((end/2)+1:end, :, :), 3).', size_data(g,1)/2, size_data(g,2), opts_expansion);

      [junk, relative_fraction{g}] = synchronize_domains(f, f);
      tmp_gfrac = f*frac_width;
      tmp_gfrac(isnan(tmp_gfrac)) = 0;
      gfraction{g} = tmp_gfrac;

      if (strncmp(fitting.scale_type, 'normalize', 10))
        %fitting.ground_truth = normalize_domain(fitting.ground_truth, f*frac_width, opts_expansion);
        fitting.ground_truth{g} = normalize_domain(fitting.ground_truth{g}, f*frac_width, opts_expansion, is_data, fitting.normalize_smooth);
        %fitting.ground_truth = normalize_domain(fitting.ground_truth, true);

        normalization_done = true;
      end

    if (strncmp(fitting.aligning_type, 'domain_width', 12))
      frac_indx(g) = find(f > fitting.fraction, 1, 'first');
    end

    if (~normalization_done & strncmp(fitting.scale_type, 'normalize', 10))
      opts_expansion = load_parameters(get_struct('ASSET'), 'domain_expansion.txt');
      [f, frac_width, full_width] = domain_expansion(mymean(fitting.ground_truth{g}(1:end/2, :, :), 3).', mymean(fitting.ground_truth{g}((end/2)+1:end, :, :), 3).', size_data(g,1)/2, size_data(g,2), opts_expansion);
      if (~fitting.fit_full)
        f(:) = 1;
      end
      %fitting.ground_truth = normalize_domain(fitting.ground_truth, f*frac_width, opts_expansion);
      fitting.ground_truth{g} = normalize_domain(fitting.ground_truth{g}, f*frac_width, opts_expansion, is_data, fitting.normalize_smooth);
      %fitting.ground_truth = normalize_domain(fitting.ground_truth, true);
      normalization_done = true;

      tmp_gfrac = f*frac_width;
      tmp_gfrac(isnan(tmp_gfrac)) = 0;
      gfraction{g} = tmp_gfrac;
    end

    if (strncmp(fitting.aligning_type, 'lsr', 3))
      sindx(g) = round(size(fitting.ground_truth{g}, 2) / 2);
      mean_ground_truth{g} = nanmean(fitting.ground_truth{g}, 3);
    end

    multi_data(g) = (size_data(g,3) > 1);
    if (multi_data(g))
      nlayers(g) = size_data(g,3);
    end

    if (~fitting.integrate_sigma)
      noise = estimate_noise(fitting.ground_truth{g});
      estim_sigma(g) = mymean(noise(:,2));
    end

    linear_truth{g} = fitting.ground_truth{g}(:);
    linear_goods{g} = isfinite(linear_truth{g});
    simul_pos = ([0:opts.nparticles-1] * opts.x_step).';
    penalty(g) = ((max(linear_truth{g}(linear_goods{g})))^2)*opts.nparticles;
    ndata(g) = length(fitting.x_pos{g});

    if (fit_temperatures)
      kB = 8.6173324e-5;
      E = 0.65;
      C2K = 273.15;
      diff_ratio = ((fitting.temperature(g)+C2K) / (opts.reaction_temperature+C2K));
      ratio = exp(-(E/kB)*((1/(fitting.temperature(g)+C2K)) - (1/(opts.reaction_temperature+C2K))));
      temp_scale{g} = [ones(1,2)*diff_ratio; ones(3,2)*ratio; ones(4,2)];
      flow_scale(g) = exp(-(E/kB)*((1/(fitting.temperature(g)+C2K)) - (1/(opts.flow_temperature+C2K))));
    else
      E = [];
    end
  end
  
  %{
  if (~fit_temperatures)
    linear_goods = linear_goods{1};
    if (strncmp(fitting.aligning_type, 'lsr', 3))
      mean_ground_truth{g} = mean_ground_truth{g};
    end
    fitting.ground_truth = fitting.ground_truth{1};
    fitting.x_pos = fitting.x_pos{1};
    fitting.t_pos = fitting.t_pos{1};
    gfraction = gfraction{1};
  end
  %}

  rescaling = orig_scaling;
  %rescaling = 10.^(floor(log10(ml_params)));
  %rescaling(ml_params == 0) = 1;
  ml_params = ml_params ./ rescaling;

  if (strncmp(fitting.type, 'simulation', 10))
    noiseless = fitting.ground_truth;
    for g=1:ngroups
      range_data(g) = fitting.data_noise * range(linear_truth{g});
    end
  end

  warning off;
  useful_data = [];

  if (fitting.cross_measure & multi_data)
    fitting.nfits = 2*fitting.nfits;
    orig_fit = fitting;
  end

  for f = 1:fitting.nfits
    if (fitting.cross_measure)
      fitting = orig_fit;
      for g=1:ngroups
        if (mod(f, 2) == 1)
          selection = randperm(size(fitting.ground_truth{g}, 3));
          nlayers(g) = round(length(selection)/2);
          fitting.ground_truth{g} = fitting.ground_truth{g}(:,:,selection(1:nlayers(g)));
        else
          fitting.ground_truth{g} = fitting.ground_truth{g}(:,:,selection(nlayers(g)+1:end));
          nlayers(g) = length(selection) - nlayers(g);
        end
        linear_truth{g} = fitting.ground_truth(:);
        linear_goods = isfinite(linear_truth{g});
      end
    end

    uuids{f} = num2str(now + cputime);
    log_name = ['adr-kymo-' uuids{f} '_'];

    tmp_params = ml_params;

    if (isempty(fitting.init_pos))
      p0 = ml_params(fit_params);
    else
      p0 = fitting.init_pos;
    end

    if (~fit_temperatures)
      flow_scale = 1;
    elseif (isempty(fitting.init_pos))
      p0 = [p0 fit_energy];
    end

    if (~fit_viscosity)
      viscosities(:) = 1;
    elseif (~all_params)
      p0 = [p0 viscosities(good_visc)];
    end

    if (fitting.fit_flow & ~all_params)
      p0 = [p0 1];
    else
      curr_flow_scale = 1;
    end

    if (fitting.fit_sigma)
      p0 = [p0 estim_sigma];
    end

    %p0 = p0 .* (1+fitting.init_noise*randn(size(p0))/sqrt(length(p0)));
    correct = 0;
    tmp_params = ml_params;
    tmp_opts = opts;

    for t=1:20
      tmp_p = real(exp(log(p0) + fitting.init_noise*randn(size(p0))/sqrt(length(p0))));
      tmp_params(fit_params) = tmp_p(1:nrates);

      tmp_opts.diffusion_params = tmp_params(1, :) .* rescaling(1, :);
      tmp_opts.reaction_params = tmp_params(2:end, :) .* rescaling(2:end, :);

      [tmp_pts, correct] = opts.init_func(tmp_opts, fitting.fit_relative, true);

      if (correct == 3)
        tmp_params = [tmp_opts.diffusion_params; tmp_pts] ./ rescaling;
        tmp_p(1:numel(fit_params)) = tmp_params(fit_params);

        break;
      elseif (fitting.start_with_best || ~isempty(fitting.init_pos))
        tmp_p = p0;
        tmp_params = ml_params;

        break;
      end
    end

    p0 = tmp_p;

    tmp_opts.diffusion_params = tmp_params(1, :) .* rescaling(1, :);
    tmp_opts.reaction_params = tmp_params(2:end, :) .* rescaling(2:end, :);

    [x0] = opts.init_func(tmp_opts, fitting.fit_relative);

    if (correct < 3)
      warning on;
      warning('Could not identify a bistable initial condition');
      warning off;
    end

    ndecimals = -min(log10(fitting.tolerance/10), 0);

    for g=1:ngroups
      if (fitting.estimate_n)
        ns = identify_n(fitting.ground_truth{g});
        nobs(g) = sum([ns(:,1); ns(:,3)]);
        norm_coeff(g) = sum(linear_goods{g})/nobs(g);
      else
        nobs(g, 1:2) = [sum(linear_goods{g}) size_data(g,2)];
        norm_coeff(g) = 1;
      end

      %fitting.score_weights = 1/nobs(1);
      curr_score_weights(g) = fitting.score_weights*length(gfraction{g})/nobs(g,1);
      half = opts.boundaries(2);

      if (fitting.integrate_sigma)
        %full_error = (0.5*nobs(1))*log(fitting.score_weights * penalty * size_data(2) * 10) + (0.5*nobs(2))*log(size_data(2)*10);
        if (fitting.pixels_only)
          full_error(g) = (0.5*nobs(g,1))*log(penalty(g) * size_data(g,2) * 10);
        else
          full_error(g) = (0.5*nobs(g,2))*log(curr_score_weights(g)*penalty(g)*size_data(g,2)*10 + size_data(g,2)*10);
        end
      else
        %full_error = (fitting.score_weights * penalty * size_data(2) * 10 + prod(size_data)) / (2*estim_sigma.^2);
        full_error(g) = (curr_score_weights(g) * penalty(g) * size_data(g,2) * 10 + prod(size_data(g,1:2))) / (2*estim_sigma(g).^2);
      end

      if (strncmp(fitting.type, 'simulation', 10))
        fitting.ground_truth{g} = noiseless{g} + range_data(g)*randn(size_data(g,:));
      end
    end
    full_error = sum(full_error);

    if (strncmp(fitting.aligning_type, 'fitting', 7) && length(p0) <= (nrates + numel(fit_energy) + fitting.fit_flow + nvisc*fit_viscosity))
      nparams = length(p0);
      fitting.aligning_type = 'domain';
      [junk, offsets] = error_function(p0(:));
      fitting.aligning_type = 'fitting';
      p0 = [p0(1:nrates) offsets/fitting.offset_scaling p0(nrates+1:end)];
    end

    nparams = length(p0);
    if (fitting.fit_flow)
      display(['Fitting ' num2str(nparams) ' parameters (' num2str(fit_params) ' & flow):']);
    else
      display(['Fitting ' num2str(nparams) ' parameters (' num2str(fit_params) '):']);
    end
    %p0 = sqrt(p0(:));

    tmp_fit = fitting;
    tmp_fit.ground_truth = [];
    tmp_fit.flow_size = size(flow);
    tmp_fit.rescale_factor = rescaling(fit_params);
    fid = fopen([log_name 'evol.dat'], 'w');
    print_all(fid, tmp_fit);
    fclose(fid);

    switch (fitting.fitting_type)
      case 'simplex'
        options = optimset('MaxFunEvals', fitting.max_iter, ...
                           'MaxIter', fitting.max_iter, ...
                           'Display', 'off', ...
                           'OutputFcn', @log_simplex, ...
                           'TolFun', fitting.tolerance, ...
                           'TolX', fitting.tolerance/10);

        uuid = ['SPLX' num2str(ceil(rand(1)*100)) ' '];
        cmaes_count = 0;

        fid = fopen([log_name 'evol.dat'], 'a');

        fprintf(fid, [uuid '%% columns="iteration, evalutation | lastbest" (' num2str(clock, '%d/%02d/%d %d:%d:%2.2f') ')\n']);

        display('Local optimization using the Simplex method');
        [p, fval_opt, stopflag_opt, out_opt] = fminsearch(@error_function, p0(:), options);
        fclose(fid);

      case 'cmaes'

        opt = cmaes('defaults');
        opt.MaxFunEvals = fitting.max_iter;
        opt.TolFun = fitting.tolerance;
        opt.TolX = fitting.tolerance/10;
        opt.SaveFilename = '';
        opt.SaveVariables = 'off';
        opt.EvalParallel = 'yes';
        opt.LogPlot = 0;
        opt.LogFilenamePrefix = log_name;
        opt.StopOnWarnings = false;
        opt.WarnOnEqualFunctionValues = false;
        opt.PopSize = 5*numel(p0);

        if (fit_temperatures)
          opt.PopSize = opt.PopSize + nvisc;
        end

        [p, fval, cmaes_count, stopflag, out] = cmaes(@error_function, p0(:), fitting.step_size, opt); 
        p = out.solutions.bestever.x;
        fval = out.solutions.bestever.f;

        options = optimset('MaxFunEvals', max(fitting.max_iter - cmaes_count,0), ...
                           'MaxIter', max(fitting.max_iter - cmaes_count,0), ...
                           'Display', 'off', ...
                           'OutputFcn', @log_simplex, ...
                           'TolFun', fitting.tolerance, ...
                           'TolX', fitting.tolerance/10);
        uuid = out.uuid;
        fid = fopen([log_name 'evol.dat'], 'a');

        display('Local optimization using the Simplex method');
        [p, fval_opt, stopflag_opt, out_opt] = fminsearch(@error_function, p(:), options);
        fclose(fid);

        disp(['Simplex: ' num2str(fval) ' -> ' num2str(fval_opt)]);

      case 'pso'
        opt = [1 2000 24 0.5 0.5 0.7 0.2 1500 fitting.tolerance 250 NaN 0 0];
        opt(2) = fitting.max_iter;
        opt(9) = 1e-10;
        opt(12) = 1;
        opt(13) = 1;
     

        bounds = [zeros(nparams, 1) ones(nparams, 1)*10];
        [pbest, tr, te] = pso_Trelea_vectorized(@error_function, nparams, NaN, bounds, 0, opt, '', p0.', [log_name 'evol']);
        p = pbest(1:end-1);
      case 'godlike'
        opt = set_options('Display', 'on', 'MaxIters', fitting.max_iter, 'TolFun', fitting.tolerance, 'LogFile', [log_name 'evol']);
        which_algos = {'DE';'GA';'PSO';'ASA'};
        which_algos = which_algos(randperm(4));
        [p, fval] = GODLIKE(@error_function, 20, zeros(nparams,1), ones(nparams,1) * 10, which_algos, opt);

      case 'dram'
%        [estim_error, estim_sigma2] = error_function(p0(:));

%        keyboard

        model.ssfun    = @error_function;

        params.par0    = p0(:); % initial parameter values
%        params.n       = nobs;
%        params.n0      = params.n;  % prior for error variance sigma^2
%        params.sigma2  = estim_sigma2;  % prior accuracy for sigma^2

        options.nsimu    = fitting.max_iter;               % size of the chain
        options.adaptint = adaptint;            % adaptation interval
        options.drscale  = drscale;

        if (numel(fitting.step_size)==nparams)
          options.qcov     = diag(fitting.step_size); % initial proposal covariance 
        elseif (numel(fitting.step_size)==(nparams^2))
          options.qcov     = fitting.step_size; % initial proposal covariance 
        else
          options.qcov     = fitting.step_size(1)*eye(nparams); % initial proposal covariance 
        end
        options.qcov = options.qcov.^2;
        %options.qcov     = fitting.step_size*eye(nparams).*2.4^2./nparams;      % initial proposal covariance 
        options.ndelays  = fitting.ndelays;
        options.stall_thresh = fitting.stall_thresh;
        options.log_file = [log_name 'evol'];
        options.printint = 10;

        % run the chain
        results = dramrun(model,[],params,options);
        p = results.mean;
      case 'sample'
        options.sampling_type = fitting.combine_data;
        options.range_size  = fitting.step_size;
        options.max_iter    = fitting.max_iter;
        options.is_log = fitting.sample_log;
        options.log_file = [log_name 'evol'];
        options.printint = 10;
        options.precision = ndecimals;

        p = exhaustive_sampler(@error_function, p0(:), options);
      otherwise
        error([fitting.fitting_type ' machine learning algorithm is not implemented']);

        return;
    end

    p = roundn(p, -(ndecimals+2));
    p(1:nrates) = abs(p(1:nrates));
    %p = p.^2;

    display(['Best (' num2str(p(:).') ')']);
  end

  warning on;

  return;

  function stop = log_simplex(x, optimValues, state)

    stop = false;

    if ((mod(optimValues.iteration, 10) == 1 && strncmp(state, 'iter', 4)) || strncmp(state, 'done', 4))
      print_str = [' %.' num2str(ndecimals) 'f'];

      if (uuid(1) == 'C')
        fprintf(fid, [uuid '1 %ld  %e 0 0'], optimValues.iteration+cmaes_count, optimValues.fval);
        fprintf(fid, print_str, x);
        fprintf(fid, '\n');
      else
        fprintf(fid, [uuid '%ld : %ld |'], optimValues.iteration+cmaes_count, optimValues.fval);
        fprintf(fid, print_str, x);
        fprintf(fid, '\n');
      end

      disp([num2str(optimValues.iteration+cmaes_count) ' : ' num2str(optimValues.fval), ' | ' num2str(x(:).')]);
    end

    return;
  end

  %function [err_all, sigma2] = error_function(varargin)
  function [err_all, offsets] = error_function(varargin)

    %sigma2 = 0;
    %p_all = varargin{1};
    p_all = roundn(varargin{1}, -ndecimals);

    [curr_nparams, nevals] = size(p_all);
    flip = false;
    if (nevals == nparams)
      flip = true;
      p_all = p_all.';
      nevals = curr_nparams; 
    end
    err_all = NaN(ngroups, nevals);

    %p_all = p_all.^2;
    %p_all = abs(p_all);

    offsets = zeros(1, ngroups);

    for i = 1:nevals
      correct = 3;

      for g=1:ngroups
        curr_p = p_all(:,i);

        tmp_params = ml_params;
        tmp_params(fit_params) = abs(curr_p(1:nrates));
        more_params = curr_p(nrates+1:end);

        if (fitting.fit_sigma)
          estim_sigma = abs(more_params(end));
          more_params = more_params(1:end-1);
        end

        if (fitting.fit_flow)
          curr_flow_scale = abs(more_params(end));
          more_params = more_params(1:end-1);
        else
          curr_flow_scale = 1;
        end

        if (fit_viscosity)
          curr_visc = ones(1,ntemps);
          curr_visc(good_visc) = abs(more_params(end-nvisc+1:end));
          more_params = more_params(1:end-nvisc);
        else
          curr_visc = viscosities;
        end

        if (fit_temperatures)
          if (isempty(fit_energy))
            tmp_params = tmp_params .* temp_scale{g};
            curr_flow_scale = curr_flow_scale * flow_scale(g);
          else
            if (fitting.fit_model)
              switch nenergy
                case 1
                  E = ones(3,2)*abs(more_params(end));
                  flow_E = E(1);
                  more_params = more_params(1:end-1);
                case 2
                  E = ones(3,2)*abs(more_params(end-1));
                  flow_E = abs(more_params(end));
                  more_params = more_params(1:end-2);
                case 4
                  E = abs(more_params(end-3:end-1));
                  E = E(:);
                  E = E(:,[1 1]);
                  flow_E = abs(more_params(end));
                  more_params = more_params(1:end-4);
                case 7
                  E = abs(more_params(end-6:end-1));
                  E = E(:);
                  E = [E(1:3,1) E(4:6,1)];
                  flow_E = abs(more_params(end));
                  more_params = more_params(1:end-7);
                otherwise
                  error(['No temperature model with ' num2str(nenergy) ' parameters has been implemented.']);
              end

              diff_ratio = ((fitting.temperature(g)+C2K) / (opts.reaction_temperature+C2K));
              ratio = exp(-(E/kB)*((1/(fitting.temperature(g)+C2K)) - (1/(opts.reaction_temperature+C2K))));
              flow_ratio = exp(-(flow_E/kB)*((1/(fitting.temperature(g)+C2K)) - (1/(opts.flow_temperature+C2K))));
            else

              diff_ratio = 1;

              curr_ratios = ones(nenergy,ntemps);
              curr_ratios(:,good_visc) = reshape(abs(more_params(end-nenergy*nvisc+1:end)), nenergy, nvisc);
              more_params = more_params(1:end-nenergy*nvisc);

              eindx = temp_indx(g);

              switch nenergy
                case 1
                  ratio = ones(3,2)*curr_ratios(1,eindx);
                  flow_ratio = ratio(1);
                case 2
                  ratio = ones(3,2)*curr_ratios(1,eindx);
                  flow_ratio = curr_ratios(2, eindx);
                case 4
                  ratio = curr_ratios(1:3, [eindx eindx]);
                  flow_ratio = curr_ratios(4, eindx);
                case 7
                  ratio = [curr_ratios(1:3, eindx) curr_ratios(4:6, eindx)];
                  flow_ratio = curr_ratios(7, eindx);
                otherwise
                  error(['No temperature scaling with ' num2str(nenergy) ' parameters has been implemented.']);
              end

              if (bad_flow(eindx))
                flow_ratio = 1;
                curr_flow_scale = 1;
              elseif (good_visc(eindx))
                curr_flow_scale = 1;
              end
            end

            tmp_params = tmp_params .* [ones(1,2)*diff_ratio; ratio; ones(4,2)];
            curr_flow_scale = curr_flow_scale * flow_ratio;
          end
        end
        if (fit_viscosity)
          tmp_params(1,:) = tmp_params(1,:)*curr_visc(temp_indx(g));
        end

        if (opts.restart_init)
          opts.diffusion_params = tmp_params(1, :) .* rescaling(1, :);
          opts.reaction_params = tmp_params(2:end, :) .* rescaling(2:end, :);
          [x0, correct] = opts.init_func(opts, fitting.fit_relative);

          if (correct == 0)
            err_all(g,i) = Inf;
            continue;
          end
        end

        normalization_done = false;
        if (fitting.fit_relative)
          [res, t] = simulate_model_real(x0, tmp_params .* rescaling, opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow * curr_flow_scale, opts.user_data, opts.max_iter);
        else
          [res, t] = simulate_model_mix(x0, tmp_params .* rescaling, opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow * curr_flow_scale, opts.user_data, opts.max_iter);
        end

       % if (~isempty(useful_data))
       %   figure;imagesc(useful_data - res);
       % end
       % useful_data = res;
        %if (t(end) < opts.tmax)
        %  err_all(i) = Inf;
        %  continue;
        %end
        
        res = res((end/2)+1:end, :);
        if (opts.nparticles ~= ndata)
          res = flipud(interp1q(simul_pos, flipud(res), fitting.x_pos{g}.'));
        end
        
        [f, fwidth] = domain_expansion(res.', size(res, 1), size(res,2), opts_expansion);
        fraction = f*fwidth;
        res = [res; res];

        switch fitting.aligning_type
          case 'best'
            if (length(t) < size_data(g,2))
              cc = normxcorr2(res, mymean(fitting.ground_truth{g}, 3)); 

              [max_cc, imax] = max(cc(size(res, 1), 1:2*size(res,2)));
              corr_offset = -(imax-size(res, 2));
            elseif (all(isfinite(res)))
              cc = normxcorr2(mymean(fitting.ground_truth{g}, 3), res); 

              [max_cc, imax] = max(cc(size_data(g,1), 1:2*size_data(g,2)));
              corr_offset = (imax-size_data(g,2));
            else
              corr_offset = 0;
            end
          case 'domain'
            if (size(res,2) <= 10)
              err_all(g,i) = Inf;
              continue;
            else

              %[f, fwidth] = domain_expansion(res(1:end/2, :).', size(res, 1)/2, size(res,2), opts_expansion);
              if (strncmp(fitting.scale_type, 'normalize', 10))
                %res = normalize_domain(res, f*fwidth, opts_expansion, false);
                res = normalize_domain(res, fraction, opts_expansion, false, fitting.normalize_smooth);
                %res = normalize_domain(res, false);
                normalization_done = true;
              end
            end

            if (isnan(f(end)))
              err_all(g,i) = Inf;
              continue;
            end
            %findx = find(f > fitting.fraction, 1, 'first');
            %corr_offset = frac_indx(g) - findx + 1;
            
            corr_offset = synchronize_domains(relative_fraction{g}, f);
            disp_args = corr_offset;
            corr_offset = corr_offset(1);

            if (isempty(corr_offset))
              corr_offset = NaN;
            end

            offsets(g) = corr_offset;
          case 'domain_width'
            if (size(res,2) <= 10)
              err_all(g,i) = Inf;
              continue;
            else

              %[f, fwidth] = domain_expansion(res(1:end/2, :).', size(res, 1)/2, size(res,2), opts_expansion);
              if (strncmp(fitting.scale_type, 'normalize', 10))
                %res = normalize_domain(res, f*fwidth, opts_expansion, false);
                res = normalize_domain(res, fraction, opts_expansion, false, fitting.normalize_smooth);
                %res = normalize_domain(res, false);
                normalization_done = true;
              end
            end

            if (isnan(f(end)))
              err_all(g,i) = Inf;
              continue;
            end
            findx = find(f > fitting.fraction, 1, 'first');
            %corr_offset = frac_indx(g) - findx + 1;
            corr_offset = findx - frac_indx(g);
          case 'fitting'
            corr_offset = more_params(g)*fitting.offset_scaling;
          case 'end'
            corr_offset = size_data(g,2) - size(res, 2) + 1;
          case 'lsr'
            rindx = min(sindx(g), size(res, 2));

            [junk, junk2, rindx] = find_min_residue(mean_ground_truth{g}, sindx(g), res, rindx, 0.95);
            corr_offset = sindx(g) - rindx;
        end

%        fraction(isnan(fraction)) = 0;
%        gindxs = [0:size(res, 2)-1];
%        goods = ((gindxs + corr_offset) > 0 & (gindxs + corr_offset) <= size_data(g,2));

%        if (sum(goods) < 10)
        %if (length(relative_fraction{g}) - corr_offset < 10)
        if (size_data(g,2) - corr_offset < 10)
          err_all(g,i) = Inf;
        else
          if (~normalization_done & strncmp(fitting.scale_type, 'normalize', 10))
            %[f, fwidth] = domain_expansion(res(1:end/2, :).', size(res, 1)/2, size(res,2), opts_expansion);
            if (~fitting.fit_full)
              fraction(:) = fwidth;
            end
            %res = normalize_domain(res, f*fwidth, opts_expansion);
            res = normalize_domain(res, fraction, opts_expansion, false, fitting.normalize_smooth);
            %res = normalize_domain(res, false);
            normalization_done = true;
          end

          fraction(isnan(fraction)) = 0;
          res(isnan(res)) = 0;
          
          res = interp2([res; fraction.'], [1:length(gfraction{g})].'+corr_offset, [1:size(fitting.ground_truth{g},1)+1]);
          fraction = res(end,:).';
          res = res(1:end-1,:);
          
          fraction(isnan(fraction)) = 0;
          res(isnan(res)) = 0;

%          indxs = corr_offset+gindxs(goods);
%          tmp_fraction = zeros(1, size_data(g,2));
%          tmp_fraction(indxs) = fraction(gindxs(goods) + 1);
          %tmp_fraction(indxs(end)+1:end) = tmp_fraction(indxs(end));

%          tmp = ones(size_data(g,1:2))*res(1, 2);
%          tmp(:, indxs) = res(:, gindxs(goods) + 1);

%          fraction = tmp_fraction(:);
%          res = tmp;
          if (multi_data)
            res = repmat(res, [1 1 nlayers(g)]);
          end

          if (~normalization_done & strncmp(fitting.scale_type, 'best', 4))
            gres = ~isnan(res(:)) & linear_goods{g};
            c = [ones(sum(gres), 1), res(gres)] \ linear_truth{g}(gres);
            
            if (c(2) <= 0)
              err_all(g,i) = Inf;
              continue;
            end

            res = c(1) + c(2)*res;
          end
          
          %tmp_err = fitting.score_weights*(fitting.ground_truth - res).^2 / norm_coeff;

          tmp_err = (fitting.ground_truth{g} - res).^2 / norm_coeff(g);
          %tmp_err = 0*tmp_err;

          %tmp_frac = ((gfraction{g} - fraction)/half).^2;
          tmp_frac = 1./(1+exp(-15*(((gfraction{g} - fraction)/half).^2-0.5)));
          %tmp_frac = 0;

          %[sum(tmp_err(linear_goods)) sum(tmp_frac)]

          %if (nargout == 2)
          %  tmp_mean = mymean(tmp_err(:));
          %  sigma2 = mymean((tmp_err(:) - tmp_mean).^2);
          %end
          %tmp_err(~linear_goods) = 0;

          if (fitting.integrate_sigma)
            % Integrated the sigma out from the gaussian error function
            %err_all(i) = (0.5*nobs(1))*log(sum(tmp_err(linear_goods))) + (0.5*nobs(2))*log(sum(tmp_frac));
            %err_all(g,i) = curr_score_weights(g)*(0.5*nobs(g,1))*log(sum(tmp_err(linear_goods{g}))) + (0.5*nobs(g,2))*log(sum(tmp_frac));

            if (fitting.pixels_only)
              err_all(g,i) = (0.5*nobs(g,1))*log(sum(tmp_err(linear_goods{g})));
            else
              err_all(g,i) = (0.5*nobs(g,2))*log(curr_score_weights(g)*sum(tmp_err(linear_goods{g})) + sum(tmp_frac));
            end
          else
            if (fitting.fit_sigma)
              %err_all(i) = sum(tmp_err(linear_goods) + sum(tmp_frac)) / (2*estim_sigma.^2) + nobs*log(estim_sigma) ;
              err_all(g,i) = (curr_score_weights(g)*sum(tmp_err(linear_goods{g})) + sum(tmp_frac)) / (2*estim_sigma(g).^2) + nobs(g)*log(estim_sigma) ;
            else
              %err_all(i) = sum(tmp_err(linear_goods) + sum(tmp_frac)) / (2*estim_sigma.^2);
              err_all(g,i) = (curr_score_weights(g)*sum(tmp_err(linear_goods{g})) + sum(tmp_frac)) / (2*estim_sigma(g).^2);
            end
          end
          % Standard log likelihood with gaussian prob
          %err_all(i) = sum(tmp_err(:) / (2*error_sigma^2));

          if (fitting.display)
            tmp_plot_avg = mymean(res,3).';
            tmp_plot_err = mymean(tmp_err,3).';
            tmp_plot_truth = mymean(fitting.ground_truth{g},3).';

            tmp_plot_frac = fraction;
            tmp_plot_frac(1:find(tmp_plot_frac, 1, 'first')-1) = NaN;
            tmp_plot_gfrac = gfraction{g};
            tmp_plot_gfrac(1:find(tmp_plot_gfrac, 1, 'first')-1) = NaN;

            yindx = find(tmp_plot_gfrac>0, 1, 'first');
            y_tick = unique([fliplr([yindx:-20:1]) yindx:20:size(tmp_plot_avg,1)]);
            y_labels = (y_tick - yindx)*10;

            xfrac = [(size(tmp_plot_err, 2)/2)-2*tmp_plot_frac(end:-1:1); (size(tmp_plot_err, 2)/2)+2*tmp_plot_frac];
            xgfrac = [(size(tmp_plot_err, 2)/2)-2*tmp_plot_gfrac(end:-1:1); (size(tmp_plot_err, 2)/2)+2*tmp_plot_gfrac];
            yfrac = [size(tmp_plot_err, 1):-1:1 1:size(tmp_plot_err,1)].';

            tmp_pos_tick = [0:50:(size(tmp_plot_err,2)/2 - 1)];
            tmp_pos_label = fitting.x_pos{1}(tmp_pos_tick + 1);
            tmp_pos_tick = [-tmp_pos_tick(end:-1:2) tmp_pos_tick] + (size(tmp_plot_err,2)/2);
            tmp_pos_label = [-tmp_pos_label(end:-1:2) tmp_pos_label];

            %figure(hfig(g));
            hax = subplot(1,3,1, 'Parent', hfig(g), 'NextPlot', 'replace');
            %hold off;
            imagesc([tmp_plot_avg(:,1:end/2) tmp_plot_avg(:,end:-1:(end/2)+1)], 'Parent', hax);
            set(hax, 'XTick', tmp_pos_tick, 'XTickLabel', tmp_pos_label, ...
                     'YTick', y_tick, 'YTickLabel', y_labels, ...
                     'NextPlot', 'add')
            colormap(hax, blueredmap)
            plot(hax, xfrac, yfrac, 'Color', [83 83 83]/255);
            title(hax, [num2str(err_all(g,i)) ' : ' num2str(sum(tmp_err(linear_goods{g}))) ', ' num2str(sum(tmp_frac))]);

            hax = subplot(1,3,2, 'Parent', hfig(g), 'NextPlot', 'replace');
            %hold off;
            %imagesc(mymean(tmp_err, 3));
            imagesc([tmp_plot_err(:,1:end/2) tmp_plot_err(:,end:-1:(end/2)+1)], 'Parent', hax);
            set(hax, 'XTick', tmp_pos_tick, 'XTickLabel', tmp_pos_label, ...
                     'YTick', y_tick, 'YTickLabel', y_labels, ...
                     'NextPlot', 'add')
            colormap(hax, redblackmap)
            %hold on;
            %plot((size(tmp_err, 2)/2)-2*gfraction{g}, 1:size(tmp_err,1), 'k');
            %plot((size(tmp_err, 2)/2)+2*gfraction{g}, 1:size(tmp_err,1), 'k');
            %plot((size(tmp_err, 2)/2)-2*fraction, 1:size(tmp_err,1), 'w');
            %plot((size(tmp_err, 2)/2)+2*fraction, 1:size(tmp_err,1), 'w');
            plot(hax, xfrac, yfrac, 'Color', [83 83 83]/255);
            plot(hax, xgfrac, yfrac, 'Color', [83 83 83]/255);

            disp_params = tmp_params .* rescaling;
            disp_params = disp_params(fit_params);

            if (fitting.fit_relative & numel(disp_params) >= 4)
              disp_params([1 3]) = disp_params([1 3]) .* disp_params([4 2]);
              %disp_params(1) = bsxfun(@times, disp_params(1), disp_params(4));
              %disp_params(3) = bsxfun(@rdivide, bsxfun(@times, disp_params(3), disp_params(2)), (1.56.^disp_params(4)));
            end
            if (numel(disp_params) > 6)
              tmp_params = NaN(ceil((numel(disp_params)+numel(E)+2)/6),6);
              tmp_params(1:numel(disp_params)) = disp_params;
              tmp_params(numel(disp_params)+1:numel(disp_params)+numel(E)) = E(:).';
              tmp_params(numel(disp_params)+numel(E)+1:numel(disp_params)+numel(E)+1) =  curr_flow_scale;
              tmp_params(numel(disp_params)+numel(E)+2:numel(disp_params)+numel(E)+2) =  corr_offset;
              title(hax, num2str(tmp_params));
            else
              if (strncmp(fitting.aligning_type, 'fitting', 7))
                title(hax, [num2str([disp_params E(:).' corr_offset])]);
              else
                title(hax, [num2str([disp_params E(:).'])]);
              end
              %title([num2str([disp_params E(:).']) ':' num2str(disp_args)]);
            end

            hax = subplot(1,3,3, 'NextPlot', 'replace', 'Parent', hfig(g));
            %hold off;
            imagesc([tmp_plot_truth(:,1:end/2) tmp_plot_truth(:,end:-1:(end/2)+1)], 'Parent', hax);
            %set(gca, 'XTick', tmp_pos_tick, 'XTickLabel', tmp_pos_label);
            %set(gca, 'YTick', y_tick, 'YTickLabel', y_labels);
            set(hax, 'XTick', tmp_pos_tick, 'XTickLabel', tmp_pos_label, ...
                     'YTick', y_tick, 'YTickLabel', y_labels, ...
                     'NextPlot', 'add')
            colormap(hax, blueredmap)
            %hold on
            plot(hax, xgfrac, yfrac, 'Color', [83 83 83]/255);
            title(hax, fitting.type)

            %keyboard
            %drawnow
          end
        end
      end
      if (fitting.display)
        set(hfig, 'Visible', 'on');
        drawnow;
        keyboard
      end
    end

    err_all(isnan(err_all)) = Inf;
    err_all = sum(err_all, 1);

    if (flip)
      err_all = err_all.';
    end

    if (all(isinf(err_all)))
      % ML crashes when all values are non-numerical
      err_all(1) = full_error*sign(err_all(1));
    end

    return;
  end
end

function domain = normalize_domain(domain, path, opts, has_noise, do_min_max)

  prct_thresh = 5;
  path = path/opts.quantification.resolution;
  [h, w, f] = size(domain);
  h = h/2;
  pos_mat = repmat([1:h].', 1, w);
  mask = bsxfun(@le, pos_mat, path.');
  mask = repmat(flipud(mask), [2, 1]);

  nplanes = size(domain, 3);

  if (do_min_max && has_noise)
    noise = estimate_noise(domain);
    domain = min_max_domain(domain, path, 3*noise(:,2));
  end

  for i=1:nplanes
    img = domain(:,:,i);
    min_val = prctile(img(~mask), prct_thresh);
    max_val = prctile(img(mask), 100-prct_thresh);
    domain(:,:,i) = (img - min_val) / (max_val - min_val);
  end

  return;
end
