function opts = find_parameters(mymovie, trackings, opts)

  params = opts.segmentation_parameters;

  p0 = gather_parameters(params, opts);
  nparams = length(p0);
  log_name = ['ML-' num2str(opts.uuid) '-' mymovie.experiment];

  dp = 1e-25;
  lbound = 0;
  ubound = 1;

  switch opts.do_ml
    case 'cmaes'
      opt = cmaes('defaults');
      opt.MaxFunEvals = 100000;
      opt.TolFun = dp;
      %opt.LBounds = zeros(nparams,1) + lbound;
      %opt.UBounds = zeros(nparams,1) + ubound;
      opt.SaveFilename = '';
      opt.LogFilenamePrefix = log_name;
      opt.EvalParallel = 'yes';

      if (strncmp(opts.parse_frames,'random',6))
        opt.Noise.on = true;
      end
      if (opts.verbosity < 3)
        opt.LogPlot = 0;
      end
      
      [p, fval, ncoutns, stopflag, out] = cmaes(@error_function, p0, 0.15, opt);
    case 'pso'
      opt = [1 2000 24 0.5 0.5 0.7 0.2 1500 dp 250 NaN 0 0];
      opt(2) = 10000;
      opt(8) = 5000;
      opt(9) = 1e-10;
      opt(12) = 1;
      opt(13) = 1;
    
      bounds = [zeros(nparams, 1) + lbound zeros(nparams, 1) + ubound];
      [pbest, tr, te] = pso_Trelea_vectorized(@error_function, nparams, NaN, bounds, 0, opt, '', p0(:).', [log_name 'evol']);
      p = pbest(1:end-1);
    case 'godlike'
      opt = set_options('Display', 'on', 'MaxIters', 10000, 'TolFun', dp, 'LogFile', [log_name 'evol']);
      which_algos = {'DE';'GA';'PSO';'ASA'};
      which_algos = which_algos(randperm(4));
      [p, fval] = GODLIKE(@error_function, 20, zeros(nparams,1) + lbound, zeros(nparams,1) + ubound, which_algos, opt);

    otherwise
      error [opts.do_ml ' machine learning algorithm is not implemented'];

      return;
  end

  opts.segmentation_parameters = insert_parameters(params, p, opts);

  return;

  %% A nested function to avoid having to pass mymovie and trackings back and 
  %% forth through the ML algorithms
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
      new_p = p_all(:, i);

      fix_round = (mod(new_p,1) == 0);
      new_p(fix_round) = new_p(fix_round) + dp;

      switch opts.segmentation_type
        case 'markers'
          params.markers = insert_parameters(params.markers, new_p, opts.ml_type);
        case 'dic'
          params.dic = insert_parameters(params.dic, new_p, opts.ml_type);
        case 'all'
          tmp = gather_parameters(params.dic, opts.ml_type);
          params.dic = insert_parameters(params.dic, new_p(1:length(tmp),1), opts.ml_type);
          params.markers = insert_parameters(params.markers, new_p(length(tmp)+1:end,1), opts.ml_type);
      end

      opts.recompute = false;
      opts.auto_save = false;
      opts.segmentation_parameters = params;
      mymovie = segment_movie(mymovie, opts);

      [mymovie, trackings] = analyze_segmentation(mymovie, trackings, opts);

      switch opts.segmentation_type
        case 'dic'
          errors = mymovie.dic.errors;
        case 'markers'
          errors = mymovie.markers.errors;
        case 'all'
          errors = [mymovie.dic.errors mymovie.markers.errors];
      end
      switch opts.ml_type
        case 'eggshell'
          errors = errors(1,:,:,:);
        case 'cortex'
          errors = errors(2,:,:,:);
      end

      errors = errors(:);
      errors = errors(~isnan(errors));

      if (length(errors) == 0)
        err = Inf;
      else
        err = mean(errors) + std(errors) + exp(sum(lbound - new_p(new_p < lbound)) + sum(new_p(new_p > ubound) - ubound)) - 1;
      end

      err_all(i) = err;
    end

    if (flip)
      err_all = err_all.';
    end

    if (all(isinf(err_all)))
      % ML crashes when all values are non-numerical
      err_all(1) = 1e5;
    end

    return;
  end
end

function p = gather_parameters(params, type)

  p = zeros(0,1);

  if (nargin == 1)
    fields = fieldnames(params);
    for i=1:length(fields)
      if (isnumeric(params.(fields{i})) && numel(params.(fields{i})) == 1 && mod(params.(fields{i}),1) ~= 0)
        p(end+1, 1) = params.(fields{i});
      end
    end
  else
    if (isstruct(type))
      switch type.segmentation_type
        case 'dic'
          p = gather_parameters(params.dic, type.ml_type);
        case 'markers'
          p = gather_parameters(params.markers, type.ml_type);
        case 'all'
          p = gather_parameters(params.dic, type.ml_type);
          p = [p; gather_parameters(params.markers, type.ml_type)];
      end
    else
      if (strncmp(type, 'eggshell', 8) || strncmp(type, 'all', 3))
        p = [p; gather_parameters(params.eggshell_params)];
        p = [p; gather_parameters(params.eggshell_weights)];
      end
      if (strncmp(type, 'cortex', 6) || strncmp(type, 'all', 3))
        p = [p; gather_parameters(params.cortex_params)];
        p = [p; gather_parameters(params.cortex_weights)];
      end
    end
  end

  return;
end

function [params, indx] = insert_parameters(params, p, type)

  if (isnumeric(type))
    indx = type;

    fields = fieldnames(params);
    for i=1:length(fields)
      if (isnumeric(params.(fields{i})) && numel(params.(fields{i})) == 1 && mod(params.(fields{i}),1) ~= 0)
        indx = indx + 1;
        params.(fields{i}) = p(indx);
      end
    end
  elseif (isstruct(type))
    switch type.segmentation_type
      case 'markers'
        params.markers = insert_parameters(params.markers, p, type.ml_type);
      case 'dic'
        params.dic = insert_parameters(params.dic, p, type.ml_type);
      case 'all'
        tmp = gather_parameters(params.dic, type.ml_type);
        params.dic = insert_parameters(params.dic, p(1:length(tmp),1), type.ml_type);
        params.markers = insert_parameters(params.markers, p(length(tmp)+1:end,1), type.ml_type);
    end
  else
    indx = 0;

    if (strncmp(type, 'eggshell', 8) || strncmp(type, 'all', 3))
      [params.eggshell_params, indx] = insert_parameters(params.eggshell_params, p, indx);
      [params.eggshell_weights, indx] = insert_parameters(params.eggshell_weights, p, indx);
    end
    if (strncmp(type, 'cortex', 6) || strncmp(type, 'all', 3))
      [params.cortex_params, indx] = insert_parameters(params.cortex_params, p, indx);
      [params.cortex_weights, indx] = insert_parameters(params.cortex_weights, p, indx);
    end
  end

  return;
end
