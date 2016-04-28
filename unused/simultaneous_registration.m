function [window, params, score] = simultaneous_registration(imgs, centers)

  nimgs = length(imgs);
  has_intuition = (nargin==2);
  if (~has_intuition || numel(centers) ~= nimgs)
    centers = ones(nimgs, 1);
  end
  
  sizes = NaN(nimgs, 3);
  for i=1:nimgs
    sizes(i,1) = size(imgs{i},1);
    sizes(i,2) = size(imgs{i},2);
    sizes(i,3) = size(imgs{i},3);
  end
  
  goods = all(sizes > 0, 2);
  sizes = sizes(goods, :);
  imgs = imgs(goods);
  centers = centers(goods);
  nimgs = length(imgs);

  window_size = min(sizes(:,1:2), [], 1);
  half = window_size(1);
  window_size(1) = window_size(1)*2;
  window = NaN([window_size sum(sizes(:,3))]);
  ntotal = numel(window);
  min_counts = 3;

  centered = imgs;
  dw = round((sizes(:,2) - window_size(2))/2);
  for i=1:nimgs
    centered{i} = imgs{i}(:,dw(i)+1:dw(i)+window_size(2),:);
  end

  min_counts = min(min_counts, sum(sizes(:,3)));

  penalty = 0;
  [err, vars] = error_function(zeros(nimgs, 1));
  %penalty = (mean(vars(isfinite(vars))) + median(vars(isfinite(vars)))) / 2;
  penalty = mean(vars(isfinite(vars))) + std(vars(isfinite(vars)));
  %[err, vars] = error_function(centers ./ sizes(:,1));

  opt = cmaes('defaults');
  opt.MaxFunEvals = 2000*nimgs;
  opt.TolFun = 1e-8;
  opt.TolX = 1e-8;
  opt.SaveFilename = '';
  opt.SaveVariables = 'off';
  opt.EvalParallel = true;
  opt.LogPlot = 0;
  opt.LogFilenamePrefix = '';
  opt.LBounds = 0;
  opt.UBounds = 1;
  opt.StopOnWarnings = false;
  opt.WarnOnEqualFunctionValues = false;
  opt.PopSize = 5*nimgs;
  opt.LogModulo = 0;

  nstart = 5;
  params = NaN(nstart, length(centers)+1);

  for i=1:nstart
    if (has_intuition)
      [params_c, fval, ncoutns, stopflag, out, bests] = cmaes(@error_function, centers(:) ./ sizes(:,1), 0.5 / i, opt); 
    else
      init = rand(length(centers), 1);
      [params_c, fval, ncoutns, stopflag, out, bests] = cmaes(@error_function, init, 0.5, opt); 
    end
    params(i,:) = [ceil(bests.x(:) .* sizes(:,1)).' bests.f];
  end
  params(params < 1) = 1;

  [score, indx] = min(params(:,end));
  params = params(indx, 1:end-1).';

  window = stack_images(imgs, params);

  return;

%  keyboard
%  window(:) = NaN;
%  for j=1:nimgs
%    window(1:min((sizes(j, 1)-params_c(j)+1), end),:,j) = centered{j}(params_c(j):min(params_c(j)+window_size(1)-1,end),:);
%  end
  %window2 = window;
  %  window2(1:min((sizes(j, 1)-params_c2(j)+1), end),:,j) = centered{j}(params_c2(j):min(params_c2(j)+window_size(1)-1,end),:);


%{
  opt = set_options('Display', 'on', 'MaxIters', 2000, 'TolFun', 1e-8, 'LogFile', '');
  which_algos = {'DE';'GA';'PSO';'ASA'};
  which_algos = which_algos(randperm(4));
  [params_g, fval] = GODLIKE(@error_function, 20, ones(size(centers)), ones(size(centers)) * sqrt(max(sizes(:,1))), which_algos, opt);
  params_g = round(params_g.^2);
  params_g(params_g<1) = 1;

  model.ssfun    = @error_function;

  params.par0    = sqrt(centers(:)); % initial parameter values
  params.n       = nimgs;

  options.nsimu    = 2000;               % size of the chain
  options.adaptint = 200;            % adaptation interval
  options.drscale  = 5;
  options.qcov     = eye(6)*0.5;      % initial proposal covariance 
  options.ndelays  = 5;
  options.stall_thresh = 0;
  options.log_file = '';

  % run the chain
  [results, chain] = dramrun(model,[],params,options);
  [junks, mins] = min(chain(:,end));
  params_d = chain(mins, 1:end-1);
  params_d = round(params_d.^2);
  params_d(params_d<1) = 1;

  params_m = results.mean;
  params_m = round(params_m.^2);
  params_m(params_m<1) = 1;

  errs(1) = error_function(sqrt(params_c));
  errs(2) = error_function(sqrt(params_g));
  errs(3) = error_function(sqrt(params_d));
  errs(4) = error_function(sqrt(params_m));

  [v,i] = min(errs)
  %}

  return;

  function [err, variance] = error_function(varargin)

    p_all = varargin{1};

    if (size(p_all,1)~=nimgs)
      p_all = p_all.';
    end
    p_all = ceil(bsxfun(@times, p_all,sizes(:,1)));
    p_all(p_all < 1) = 1;

    for n=1:size(p_all, 2)
      p = p_all(:,n);

      window(:) = NaN;
      count = 1;
      for j=1:nimgs
        window(max(half-p(j)+2, 1):min(half+(sizes(j, 1)-p(j)+1), end),:,count:count+sizes(j,3)-1) = centered{j}(max(p(j)-half, 1):min(p(j)+half-1,end),:,:);
        count = count + sizes(j,3);
      end

      %avg = mymean(window, 3);
      avg = nanmean(window, 3);
      counts = sum(~isnan(window), 3);
      avg(counts < min_counts) = NaN;

      variance = bsxfun(@minus, window, avg).^2;
      goods = ~isnan(variance);
      %err = -(sum(variance(goods)) + log(sum(goods(:))/ntotal));
      err(n) = sum(variance(goods)) + penalty*(sum(~goods(:)));

      %err(n) = err(n)*mean(1 + p(:)./sizes(:,1));
      %[sum(variance(goods)) penalty*(sum(~goods(:))) mean(1 + p(:)./sizes(:,1)) err(n)]

      if (~any(goods))
        err(n) = NaN;
      end
    end

    return;
  end
end
