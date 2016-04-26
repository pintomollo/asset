function study_kernel_radius(radius, niter, nfirst)

  scaling = 1;
  scaling_std = 1;

  %nparticles = [30 40 50:50:800];
  nparticles = [30:10:800];
  if (nargin == 0)
    radius = [0.1:0.1:1 2 4 0.8333 0.895];
    niter = 5000;
    nfirst = 15;
  end
  aim_indx = find(nparticles>450, 1, 'first') - 1;

  npart = length(nparticles);
  nradius = length(radius);

  colors = redbluemap(nradius+1);
  colors = colors([3:end-1 1 end 2],:);

  mkernel = NaN(npart, nradius);
  ekernel = NaN(npart, nradius);
  skernel = NaN(npart, nradius);

  mavg = NaN(npart, 1);
  savg = NaN(npart, 1);

  nx = [0:2:14];

  data = load('niwayama.mat');
  fittedmodel = data.fittedmodel;
  fittedstds = data.fittedstds;

  goal = feval(fittedmodel, nx(1))*scaling;

  %nparticles = 2000;
  for n = 1:npart
    prob_size = [nparticles(n), niter];
    xr = rand(prob_size)*range(nx)+nx(1);
    yr = randn(prob_size) .* reshape(feval(fittedstds, xr)*scaling_std, prob_size) + reshape(feval(fittedmodel, xr)*scaling, prob_size);

    for r = 1:nradius
      weight = @(x)(exp(-(x.^2)/(2*(radius(r)/2)^2)));

      weights = weight(xr);
      weights = bsxfun(@rdivide, weights, sum(weights, 1));

      kernel = (sum(yr .* weights, 1) - goal);

      ekernel(n,r) = mymean(kernel);
      [mkernel(n,r), skernel(n,r)] = mymean(abs(kernel));
    end

    [xr, indx] = sort(xr, 1);
    for j = 1:prob_size(2)
      yr(:,j) = yr(indx(:,j),j);
    end

    avg = abs(mean(yr(1:min(end, nfirst), :), 1) - goal);
    [mavg(n,1), savg(n,1)] = mymean(avg);

    disp([num2str(n) '/' num2str(npart)]);
  end

  legends = cellstr([repmat('kernel ', nradius, 1) num2str(radius(:))]);
  legends{end+1} = ['average ' num2str(nfirst)];

  figure;
  hsub = subplot(2,2,1);hold on;
  set(hsub, 'ColorOrder', colors);
  plot(nparticles, mkernel(:,1:end-2))
  plot(nparticles, mkernel(:,end-1), 'LineWidth', 2, 'Color', colors(end-2,:))
  plot(nparticles, mkernel(:,end), 'LineWidth', 2, 'Color', colors(end-1,:))
  plot(nparticles, mavg, 'LineWidth', 2, 'Color', colors(end,:));
  legend(legends);
  xlabel('Particle number')
  ylabel('Average error (um/s)')

  hsub = subplot(2,2,2);hold on;
  set(hsub, 'ColorOrder', colors);
  plot(nparticles, skernel(:,1:end-2))
  plot(nparticles, skernel(:,end-1), 'LineWidth', 2, 'Color', colors(end-2,:))
  plot(nparticles, skernel(:,end), 'LineWidth', 2, 'Color', colors(end-1,:))
  plot(nparticles, savg, 'LineWidth', 2, 'Color', colors(end,:));
  xlabel('Particle number')
  ylabel('Std of error (um/s)')
  legend(legends);

  hsub = subplot(2,2,3);hold on;
  set(hsub, 'ColorOrder', colors);
  plot(nparticles, mkernel(:,1:end-2) + skernel(:,1:end-2))
  plot(nparticles, mkernel(:,end-1) + skernel(:,end-1), 'LineWidth', 2, 'Color', colors(end-2,:))
  plot(nparticles, mkernel(:,end) + skernel(:,end), 'LineWidth', 2, 'Color', colors(end-1,:))
  plot(nparticles, mavg + savg, 'LineWidth', 2, 'Color', colors(end,:));
  xlabel('Particle number')
  ylabel('Sum of both (um/s)')
  legend(legends);

  ekernel = ekernel/abs(goal);
  hsub = subplot(2,2,4);hold on;
  set(hsub, 'ColorOrder', colors);
  plot(nparticles, ekernel(:,1:end-2))
  plot(nparticles, ekernel(:,end-1), 'LineWidth', 2, 'Color', colors(end-2,:))
  plot(nparticles, ekernel(:,end), 'LineWidth', 2, 'Color', colors(end-1,:))
  legend(legends);
  title(num2str(mean(ekernel(1:aim_indx, end))));
  xlabel('Particle number')
  ylabel('Average value (%)')

  keyboard

  nparticles = [10:20:450];
  npart = length(nparticles);

  opt = cmaes('defaults');
  opt.MaxFunEvals = 20000;
  opt.TolFun = 1e-8;
  opt.TolX = 1e-8;
  opt.SaveFilename = '';
  opt.SaveVariables = 'off';
  opt.EvalParallel = true;
  opt.LogPlot = 0;
  opt.LogFilenamePrefix = '';
  opt.LBounds = 0.01;
  opt.UBounds = 1;
  opt.StopOnWarnings = false;
  opt.WarnOnEqualFunctionValues = false;
  opt.PopSize = 15;
  opt.LogModulo = 0;

  ntests = 5;
  all_bests = NaN(ntests, 2);
  for i=1:ntests
    [params_c, fval, ncoutns, stopflag, out, bests] = cmaes(@error_function, 0.75, 0.2, opt);
    all_bests(i,:) = [bests.x bests.f];
  end

  keyboard

  return;

  function [err] = error_function(varargin)

    p_all = varargin{1};
    nradius = length(p_all);
    mkernel = NaN(npart, nradius);
    skernel = NaN(npart, nradius);

    for n = 1:npart
      prob_size = [nparticles(n), niter];
      xr = rand(prob_size)*range(nx)+nx(1);
      yr = randn(prob_size) .* reshape(feval(fittedstds, xr), prob_size) + reshape(feval(fittedmodel, xr), prob_size);

      for r = 1:nradius
        weight = @(x)(exp(-(x.^2)/(2*(p_all(r)/2)^2)));

        weights = weight(xr);
        weights = bsxfun(@rdivide, weights, sum(weights, 1));

        kernel = abs(sum(yr .* weights, 1) - goal);

        [mkernel(n,r), skernel(n,r)] = mymean(kernel);
      end
    end

    err = sum(mkernel + skernel, 1);
    err = reshape(err, size(p_all));

    return;
  end
end
