function opts = find_domain(opts, uuid)

  if (nargin == 0)
    opts = get_struct('ASSET');
    opts.do_ml = 'cmaes';
  end
  %opts = load_parameters(opts, 'domains');
  opts = load_parameters(opts, 'domain_center');

  if (opts.uuid == 0)
    opts.uuid = now + cputime;
  end

  manuals = dir('*_DP.mat');
  manuals = manuals(randperm(length(manuals)));

  p0 = gather_parameters(opts.quantification, opts);

  nparams = length(p0);
  log_name = ['ML-' num2str(opts.uuid) '-domain_'];

  niter = 1;
  counts = niter+1;
  orig = [];

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
      
      [p, fval, ncoutns, stopflag, out] = cmaes(@error_function, p0, 0.35, opt);
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

  opts.quantification = insert_parameters(opts.quantification, p, opts);

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

    %if (counts > niter)
    %  counts = 0;
    %  manuals = manuals([2:end 1]);
    %  orig = load(manuals(1).name);

    %  while (isempty(orig.path))
    %    manuals = manuals([2:end]);
    %    orig = load(manuals(1).name);
    %  end
      %tmp = double(orig.domain);
      %tmp(orig.domain == intmax(class(orig.domain))) = NaN;
      %orig.domain = imnorm(tmp);
    %end

    %estim_mid = size(orig.domain, 2) / 2;
    %keyboard

    for i = 1:nevals
      new_p = p_all(:, i);
      real_p = new_p;
      real_p(real_p <= 0) = dp;
      real_p(real_p > 1) = 1 - dp;

      %fix_round = (mod(new_p,1) == 0);
      %new_p(fix_round) = new_p(fix_round) + dp;

      opts.quantification = insert_parameters(opts.quantification, real_p, opts);
      
      for m=1:length(manuals)
        
        orig = load(manuals(m).name);
        if (size(orig.path, 1) ~= size(orig.domain, 1))
          continue;
        end
        estim_mid = size(orig.domain, 2) / 2;

        path = dynamic_programming(orig.domain, opts.quantification.params, @weight_symmetry, opts.quantification.weights, opts);

        wimg = weight_symmetry(orig.domain, opts.quantification.weights);
        values = bilinear_mex(wimg, path, 1:length(path));
        indx = find(values < opts.quantification.weights.gamma, 1, 'first');

        if (isempty(indx))
          indx = length(path);
        end

        path(1:indx-1) = NaN;

        %path = dynamic_prog_2d(orig.domain, opts.quantification.params, @weight_domain, opts.quantification.weights, @init_domain, opts.quantification.init_params, opts);

        errors = abs(path(:) - orig.path(:,1));
        %errors(:, 1) = errors(:, 1) * 2;
        %errors = sum(errors, 2);
        errors(isnan(errors)) = 2*abs(path(isnan(errors), 1) - estim_mid);
        errors(isnan(errors)) = 2*abs(orig.path(isnan(errors), 1) - estim_mid);
        errors = errors(~isnan(errors));

        if (length(errors) == 0)
          err = Inf;
        else
          err = mean(errors(:)) + std(errors(:)) + exp(3*(sum(lbound - new_p(new_p < lbound)) + sum(new_p(new_p > ubound) - ubound))) - 1;
        end

        if (isnan(err_all(i)))
        err_all(i) = err;
        else
        err_all(i) = err_all(i) + err;
        end
      end
    end

    if (flip)
      err_all = err_all.';
    end

    if (all(isinf(err_all)))
      % ML crashes when all values are non-numerical
      err_all(1) = 1e5;
    end

    counts = counts + nevals;

    return;
  end
end

function p = gather_parameters(params, type)

  p = zeros(0,1);

  if (nargin == 1)
    fields = fieldnames(params);
    for i=1:length(fields)
      for j=1:length(params.(fields{i}))
        if (isnumeric(params.(fields{i})(j)) && mod(params.(fields{i})(j),1) ~= 0)
          p(end+1, 1) = params.(fields{i})(j);
        end
      end
    end
  else
    p = gather_parameters(params.init_params);
    p = [p; gather_parameters(params.params)];
    p = [p; gather_parameters(params.weights)];
  end

  return;
end

function [params, indx] = insert_parameters(params, p, type)

  if (isnumeric(type))
    indx = type;

    fields = fieldnames(params);
    for i=1:length(fields)
      for j=1:length(params.(fields{i}))
        if (isnumeric(params.(fields{i})(j)) && mod(params.(fields{i})(j),1) ~= 0)
          indx = indx + 1;
          params.(fields{i})(j) = p(indx);
        end
      end
    end
  else
    indx = 0;

    [params.init_params, indx] = insert_parameters(params.init_params, p, indx);
    [params.params, indx] = insert_parameters(params.params, p, indx);
    [params.weights, indx] = insert_parameters(params.weights, p, indx);
  end

  return;
end
