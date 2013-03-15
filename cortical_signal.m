function mymovie = cortical_signal(mymovie, opts)

  type = opts.segmentation_type;
  if (strncmp(type, 'markers', 7) | strncmp(type, 'all', 3))
    [nframes imgsize] = size_data(mymovie.cortex);
  else
    [nframes imgsize] = size_data(mymovie.(type));
  end

  if (~opts.recompute && isfield(mymovie.data, 'quantification') && ...
     ~isempty(mymovie.data.quantification) && length(mymovie.data.quantification) == nframes)

    return;
  end

  window_shape = opts.quantification.window_shape;

  if (ischar(window_shape))
    window_size = opts.quantification.window_size / opts.pixel_size;

    switch window_shape
      case 'gaussian'
        window_shape = fspecial(window_shape, ceil(window_size), opts.quantification.window_params / opts.pixel_size);

      otherwise
        window_shape = fspecial(window_shape, ceil(window_size));
    end

    window_shape = window_shape / sum(window_shape(:));
  end

  if (opts.quantification.use_ruffles)
    if ((opts.recompute & opts.segment) | ~isfield(mymovie.(type), 'ruffles') | isempty(mymovie.(type).ruffles))
    %if (~isfield(mymovie.(type), 'ruffles') | isempty(mymovie.(type).ruffles))
      mymovie = find_ruffles(mymovie, opts);
      mymovie = track_ruffles(mymovie, opts);
      mymovie = follow_invaginations(mymovie, opts);
    elseif (~isfield(mymovie.(type).ruffles, 'paths') | empty_struct(mymovie.(type).ruffles, 'paths'))
    %elseif (~isfield(mymovie.(type).ruffles, 'paths') | isempty(mymovie.(type).ruffles(1).paths))
      %do_it = true;
      %for j=1:nframes
      %  if (~isempty(mymovie.(type).ruffles(j).paths))
      %    do_it = false;
%
      %    break;
      %  end
      %end
      %if (do_it)
        mymovie = follow_invaginations(mymovie, opts);
      %end
    end
  end

  if ((opts.recompute & opts.segment) | ~isfield(mymovie.data, 'eggshell') | isempty(mymovie.data.eggshell))
    mymovie = duplicate_segmentation(mymovie, 'data', opts);
    mymovie.data = smooth_segmentation(mymovie.data, opts);
  end

  % Progress bar
  if (opts.verbosity > 0)
    hwait = waitbar(0,'Quantifying signal','Name','ASSET');
  end

  optims = optimset('Display', 'off');
  bin_dist = 64;
  bin_step = 1;
  prev_width = NaN;
  ninit = 10;

  priors = NaN(ninit, 2);
%  egg_priors = NaN(ninit, 2);
  for i=1:ninit
    nimg = randi(nframes);
    ph = double(load_data(mymovie.cortex, nimg));
%    data = double(load_data(mymovie.data, nimg));
    cortex = mymovie.data.cortex(nimg).carth;

    if (~isempty(cortex))
      [priors(i,:)] = estimate_peak(ph, cortex, [], true, opts);
%      [egg_priors(i,:)] = estimate_peak(data, mymovie.data.eggshell(nimg).carth, [0 priors(i,2)], true, opts);
    end
  end
  
  priors = median(priors, 1);
%  egg_priors = median(egg_priors, 1);
  egg_stats = NaN(nframes, 2);

  for i=1:nframes
    %bounds = [-Inf bin_step 0 0 -Inf 0 0; ...
    %           Inf Inf Inf Inf Inf Inf Inf];

    %bounds_invag = [-Inf bin_step -Inf -Inf 0 0; ...
    %                 Inf Inf Inf Inf Inf Inf];

    nimg = i;
    %nimg = randi(nframes)
    %nimg = i + 30
    %nimg = 30

    cortex = mymovie.data.cortex(nimg).carth;

    if (isempty(cortex))
      continue;
    end

    if (opts.quantification.use_ruffles)
      [ cortex, rescale] = insert_ruffles(cortex, mymovie.data.ruffles(nimg).paths);
    else
      rescale = false(size(cortex, 1), 1);
    end

    if (opts.recompute||(~isfield(mymovie.data, 'quantification'))||(length(mymovie.data.quantification) < nimg)||~isfield(mymovie.data.quantification(nimg), 'front')||isempty(mymovie.data.quantification(nimg).front))
      img = double(load_data(mymovie.data, nimg));
      img = mask_neighbors(img, mymovie.data.centers(:,nimg), mymovie.data.axes_length(:,nimg), mymovie.data.orientations(1,nimg), mymovie.data.neighbors(nimg), opts);

      ph = double(load_data(mymovie.cortex, nimg));
      ph = mask_neighbors(ph, mymovie.data.centers(:,nimg), mymovie.data.axes_length(:,nimg), mymovie.data.orientations(1,nimg), mymovie.data.neighbors(nimg), opts);

      % Background substraction

      safety = ceil(2 / opts.pixel_size);
      path = draw_ellipse(mymovie.data.centers(:, nimg), mymovie.data.axes_length(:, nimg)+safety, mymovie.data.orientations(nimg));
      mask = roipoly(img, path(:,1), path(:,2));
      img2 = img;
      img2(mask) = NaN;
      vertical_bkg = nanmean(img2,2);
      img = bsxfun(@minus, img, vertical_bkg);
      img2 = img;
      img2(mask) = NaN;
      horizontal_bkg = nanmean(img2,1);
      img = bsxfun(@minus, img, horizontal_bkg);
      img(img < 0) = 0;

      egg = mymovie.data.eggshell(nimg).carth;
      egg_dist = min(sqrt(bsxfun(@minus, egg(:,1), cortex(:,1).').^2 + bsxfun(@minus, egg(:,2), cortex(:,2).').^2), [], 2);
      %egg = mymovie.markers.eggshell(nimg).carth;
      %if (nimg>1)
      %  if (prev_width > 0)
      %    prev_width = 0.9*prev_width + 0.1*mean(mymovie.data.quantification(nimg).front(:,3));
      %  else
      %    prev_width = mean(mymovie.data.quantification(nimg).front(:,3));
      %  end
      %end

%      keyboard

%      [curr_prior, dperp, dpos] = estimate_peak(img, cortex, priors, opts);
%      curr_prior

%      if (isfinite(prev_width))
%        prev_width = 0.9*prev_width + 0.1*peak_prior(2);
%      else
%        prev_width = peak_prior(2);
%      end
      %peak_prior = priors;
      %peak_prior(1) = 0.9*priors(1) + 0.1*curr_prior(1)

      %[gs] = estimate_signal(img, cortex, dperp, dpos, peak_prior, opts);
      pos = estimate_peak(ph, cortex, priors, false, opts);
      imf = emdc([], pos(:,1), true, 4);
      %if (size(imf, 1) > 4)
      %  pos = sum(imf(end-3:end, :), 1).';
      %else
      %  pos = sum(imf(end-2:end, :), 1).';
      %end
      pos = imf(end, :).';

      goods = (egg_dist >= 3*priors(2));
      [gs, dperp] = estimate_signal(img, cortex, pos, priors(2), opts);

      if (any(goods))
        egg_pos = estimate_peak(img, egg, [0 priors(2)/2], false, opts);
      
        %%%%%%%%%%%% USEFEUL TO TRY MEAN AS WELL
        %egg_width = nanmean(egg_pos(goods, 2));
        egg_width = priors(2);

        indxs = [1:size(egg_pos, 1)].';

        imf = emdc(indxs(goods), egg_pos(goods,1), true, 4);
        egg_pos = NaN(size(indxs));
        egg_pos(goods) = imf(end, :).';

        %[egg_gs] = estimate_signal(img, egg, egg_pos, priors(2), opts);
        [egg_gs] = estimate_signal(img, egg, egg_pos, egg_width, opts);
        egg_stats(nimg, :) = nanmean(egg_gs(:, [1 3]));
      end

%      [gs, dperp] = new_estimate(img, ph, cortex, egg, prev_width, opts);

      %keyboard

      %{
      [values, dperp, dpos] = perpendicular_sampling(img, cortex, opts);
      [ph_values] = perpendicular_sampling(ph, cortex, dperp, dpos, opts);

      egg_vect_x = bsxfun(@minus, egg(:,1).', cortex(:,1));
      egg_vect_y = bsxfun(@minus, egg(:,2).', cortex(:,2));

      dotprod = bsxfun(@times, dperp(:,1), egg_vect_x) + bsxfun(@times, dperp(:,2), egg_vect_y);
      norm = (egg_vect_x.^2 + egg_vect_y.^2);
      dist = (norm - dotprod.^2);
      dist(dotprod < 0) = Inf;
      [junk, conn_indx] = min(dist, [], 2);
      egg_dist = dotprod(sub2ind(size(dotprod), [1:size(cortex, 1)].', conn_indx));

%      figure;subplot(1,2,1)
%      imagesc(values)
%      subplot(1,2,2)
%      imagesc(ph_values)
%      mtit(num2str(nimg))
%      drawnow

      npts = size(values, 1);
      slopes = NaN(npts, 2);

      [params, bounds] = estimate_mean(dpos, ph_values(~rescale, :), egg_dist(~rescale), false, bounds, optims);
      params(1:2)
      [params, junk, slopes(~rescale, 1)] = estimate_mean(dpos, values(~rescale, :), egg_dist(~rescale), false, bounds, optims);
      params(1:2)

      egg_values = perpendicular_sampling(img, egg(conn_indx, :), dperp, dpos(dpos > params(2)*1.79), opts);

      keyboard

      bounds_invag(:, 1) = bounds(:, 1) + bounds(:, 2);
      bounds_invag(:, 2) = 2*bounds(:,2);
      [params_invag, bounds_invag, slopes(rescale, :)] = estimate_mean(dpos, values(rescale, :), NaN, true, bounds_invag, optims);

      signal = NaN(npts, 7); 

      for j=1:npts

        if (rescale(j))
          curr_bounds = bounds_invag;
          curr_params = [params_invag(1:2) slopes(j, :) params_invag(end-1:end)];
        else
          curr_bounds = bounds;
          curr_params = [params(1:4) slopes(j, 1) params(end-1:end)];
        end

        line = values(j,:);
        valids = ~isnan(line);

        if (~any(valids))
          warning(['No valid data point for quantification in frame ' num2str(nimg) ', line ' num2str(j) '.'])

          continue;
        end

        valid_line = line(valids);
        valid_dpos = dpos(valids);

        tmp_val = emdc(valid_dpos, valid_line);

        if (size(tmp_val, 1) > 2)
          smoothed = sum(tmp_val(end-1:end, :), 1);
        else
          smoothed = tmp_val(end, :);
        end

        %keyboard

        line_params = estimate_front(valid_dpos, smoothed, curr_params, rescale(j));
        line_params = max(line_params, curr_bounds(1,:));
        line_params = min(line_params, curr_bounds(2,:));
        line_params = fit_front(@front, line_params, valid_dpos.', valid_line.', curr_bounds(1,:), curr_bounds(2, :), optims);

        if (rescale(j))
          line_params = [line_params(1:2) NaN line_params(3:end)];
        end

        signal(j, :) = line_params;
      end
      %}

      %mymovie.data.quantification(nimg).front = signal(:, [1 2 end]);
      %mymovie.data.quantification(nimg).bkg = signal(:, [3:end-1]);
      %mymovie.data.quantification(nimg).carth = cortex;
      mymovie.data.quantification(nimg).front = gs;
      mymovie.data.quantification(nimg).bkg = {vertical_bkg, horizontal_bkg};
      mymovie.data.quantification(nimg).ruffles = rescale;
      mymovie.data.quantification(nimg).carth = cortex;
    else

      [dperp, dpos] = perpendicular_sampling(cortex, opts);
    end

    if (opts.recompute|length(mymovie.data.quantification) < nimg|~isfield(mymovie.data.quantification(nimg), 'cortex')|isempty(mymovie.data.quantification(nimg).cortex))
      mymovie.data.quantification(nimg).cortex = extract_ridge(mymovie.data.quantification(nimg).front, cortex, dperp, rescale, opts);
    end

    if (opts.verbosity > 0)
      waitbar(i/nframes,hwait);
    end
  end

  egg_stats = median(egg_stats(~isnan(egg_stats(:,1)), :));

  for i=1:nframes
    mymovie.data.quantification(i).eggshell = egg_stats;
  end

  if (opts.recompute|~isfield(mymovie.data, 'domain')|isempty(mymovie.data.domain))
    mymovie = carth2normalized(mymovie, opts);

%    params = opts.quantification.params;
%    weights = opts.quantification.weights;
%    init_params = opts.quantification.init_params;

    %img = gather_quantification(mymovie, opts);
    %img = imnorm(img);
    %path = dynamic_prog_2d(img, params, @weight_domain, weights, @init_domain, init_params, opts);

    %mymovie.data.domain = path / size(img, 2);
  end

  if (opts.verbosity > 0)
    close(hwait);
  end

  return;
end

function [sigma, ks] = normality_test(x, y, factor)
% Kolmogorov-Smirnov test

  if (nargin == 2)
    factor = 5;
  end

  min_length = 2;
  half = find(x==0);
  max_length = floor(min(half-1, length(x)-half)/(factor)) - 1;
  ks = NaN(max_length, 1);

  if (max_length < min_length)
    sigma = NaN;
    ks = NaN;
    return;
  end

  %full_cdf = cumtrapz(y);

  for i=min_length:max_length
    %figure;
    w = x(half+i);
    pos = (abs(x) <= factor*w);

    %first = find(pos, 1, 'first')-1;
    %last = find(pos, 1, 'last')+1;
    %b = y(first);
    %slope = 0.5*(y(last)-b) / ((factor*w+1));
    %npts = sum(pos);
    %p = [0:npts+1];

    bkg = interp1(x(~pos), y(~pos), x(pos), 'linear');

    tmp_x = x(pos) / w;
    tmp_y = y(pos) - bkg;
    %tmp_y2 = y(first:last) - (slope*p + b);

    tmp_y = cumsum(tmp_y);
    tmp_y = tmp_y / tmp_y(end);
    %tmp_y2 = cumtrapz(tmp_y2);
    %tmp_y2 = tmp_y2(2:end-1) / tmp_y2(end);

    %tmp_y3 = full_cdf(pos) - full_cdf(first) - (0.5*slope*p(2:end-1).^2 + b*p(2:end-1));
    %tmp_y3 = tmp_y3/(full_cdf(last) - ((0.5*slope*(last-first).^2 + b*(last-first))) - full_cdf(first));

    norm_vals = cdf('norm', tmp_x, 0, 1);
    delta = abs(norm_vals - tmp_y);
    ks(i, 1) = max(delta);
    %delta = abs(norm_vals - tmp_y2);
    %ks(i, 2) = max(delta);
    %delta = abs(norm_vals - tmp_y3);
    %ks(i, 3) = max(delta);

    %plot(tmp_x, norm_vals, 'k');
    %hold on;
    %plot(tmp_x, tmp_y, 'b');
    %plot(tmp_x, tmp_y2, 'r');
    %plot(tmp_x, tmp_y3, 'g');
  end

  [ks, indx] = min(ks, [], 1);
  sigma = x(half+indx);

  return;
end

function [peak_params, dperp, dpos] = estimate_peak(ph, cortex, priors, do_avg, opts)

%  if (nargin == 3)
%    opts = priors;
%    priors = [];
  if (isempty(priors))
    ph = gaussian_mex(ph, 1.5);
    [ph_values, dperp, dpos] = perpendicular_sampling(ph, cortex, opts);
  else
    ph = gaussian_mex(ph, priors(2)/2);
    [ph_values, dperp, dpos] = perpendicular_sampling(ph, cortex, opts);
  end

  if (do_avg)
    values = nanmean(ph_values, 1);
  else
    values = ph_values;
  end

  nrows = size(values, 1);
  peak_params = NaN(nrows, 2);

  for n=1:nrows
    y = values(n, :);
    x = dpos;

    valids = isfinite(y);

    y = y(valids);
    x = x(valids);

    %dy = differentiator(x, y, true);
    dy = differentiator_mex(x(:), y(:), true);

    %maxs = find(dy(1:end-1) > 0 & dy(2:end) <= 0);
    maxs = x(dy(1:end-1) > 0 & dy(2:end) <= 0);

    if (~isempty(priors))
      %maxs = maxs(abs(x(maxs) - priors(1)) < 2*priors(2));
      maxs = maxs(abs(maxs - priors(1)) < 2*priors(2));
    end

    if (isempty(maxs))
      if (~isempty(priors))
        peak_params(n, :) = priors;
      end
    else
%      ks = NaN(length(maxs), 2);

%      for i=1:length(maxs)
%        [ks(i,1), ks(i,2)] = normality_test(x - x(maxs(i)), y);
%      end

%      [junk, indx] = min(ks(:,2));
%      peak_params(n, :) = [x(maxs(indx)) ks(indx, 1)];

      peak_params(n, :) = normality_test_mex(x, y, maxs);

      %hold off;
      %plot(x,y);
      %hold on;
      %scatter(x(x==peak_params(n,1)),y(x==peak_params(n,1)), 'r');
      %keyboard
    end
  end

  return;
end

function [gaussians, dperp] = estimate_signal(img, cortex, pos, peak_width, opts)
 
  [values, dperp, dpos] = perpendicular_sampling(img, cortex, opts);
  
  halfed = gaussian_mex(img, peak_width/2);
  [halfed] = perpendicular_sampling(halfed, cortex, dperp, dpos, opts);

  goods = ~isnan(pos);

  coef = 5;
  peak_dist = coef*peak_width;

  gauss3 = (abs(bsxfun(@minus, dpos, pos)) < peak_dist);
  gauss3 = (nanmean(gauss3, 1) > 0.5*mean(goods));
  gauss3([1 end]) = false;

  halfed(~goods, :) = NaN;

  [X,Y] = meshgrid(dpos, [1:size(values,1)].');

  bkg = halfed;
  bkg(:,gauss3) = interp2(X(:,~gauss3), Y(:,~gauss3), bkg(:,~gauss3), X(:,gauss3), Y(:, gauss3), 'linear');

  s = halfed - bkg;
  [maxs, indxs] = max(s(:, gauss3), [], 2);
  avgs = mymean(s(:, gauss3), 2);
  % Average over [-3 3]*sigma of a gaussian is ~0.42*ampl
  ampls = (maxs-avgs)/0.58;
  goods = ~isnan(ampls);

  if (~all(goods))
    tmp_pos = [1:length(ampls)].';
    [junk, ampls] = interp_elliptic(tmp_pos(goods), ampls(goods), tmp_pos, tmp_pos([1 end]));
  end

  imf = emdc([], ampls, true, 3);
  %if (size(imf, 1) > 4)
  %  ampls = sum(imf(end-3:end, :), 1);
  %else
    ampls = imf(end, :);
  %end

  %bkg = halfed;
  %bkg(:,gauss3) = interp2(X(:,~gauss3), Y(:,~gauss3), bkg(:,~gauss3), X(:,gauss3), Y(:, gauss3), 'linear');

  ampls(ampls < 0) = 0;
  %approx3 = bkg + bsxfun(@times, ampls.', g);

  npts = length(ampls);

  tmp = ones(npts, 1);
  p0 = [ampls(:) pos peak_width*tmp];

  lbound = [0*tmp -2*peak_width+pos peak_width*0.95*tmp];
  ubound = [max(halfed(:))*tmp 2*peak_width+pos peak_width*1.05*tmp];

  gaussians = perform_fit(p0, lbound, ubound, dpos, values, halfed, coef);

  if (opts.verbosity == 3)
%    figure;
%    imshow(realign(img,[388 591],centers(:,nimg),orientations(1,nimg)));
%    hold on;
%    myplot(realign(carths,[388 591],centers(:,nimg),orientations(1,nimg)),'g');

    half = round(size(values, 1)/2);
    figure;imagesc(values([half+1:end, 1:half], :));

    p = bsxfun(@minus, dpos, gaussians(:,2));
    signal = bsxfun(@times, gaussians(:,1), exp(bsxfun(@rdivide, -p.^2, 2*gaussians(:,3).^2)));
    figure;imagesc(signal([half+1:end, 1:half], :));
  end

  %p = bsxfun(@minus, dpos, gaussians(:,2));
  %gauss = logical(median(double(bsxfun(@lt, abs(p), gaussians(:,3)*coef)), 1));

  %bkg = halfed;
  %bkg(:,gauss) = interp2(X(:,~gauss), Y(:,~gauss), bkg(:,~gauss), X(:,gauss), Y(:,gauss), 'linear');
  %final3 = bkg + bsxfun(@times, gaussians(:,1), exp(bsxfun(@rdivide, -p.^2, 2*gaussians(:,3).^2)));

  return;
end

function gaussians = perform_fit(p0, lbound, ubound, dpos, values, smoothed, coef)
%  [maxs, indxs] = max(smoothed(:, gauss), [], 2);
%  avgs = mymean(smoothed(:, gauss), 2);
  % Average over [-1.79 1.79]*sigma of a gaussian is ~0.65*ampl
%  ampls = (maxs-avgs)/0.35;

  [niter, params] = fit_cminpack_mex(p0.', values.', coef, lbound.', ubound.', dpos, smoothed.');
  params = params.';

  goods = (params(:,3) ~= -1);
  if (~all(goods))
    tmp_pos = [1:length(goods)].';
    [junk, params] = interp_elliptic(tmp_pos(goods), params(goods, :), tmp_pos, tmp_pos([1 end]));
  end


  %%%%%%%%%%%%%%% Need to handle NaNs where params(:,3) < 0

%  bound = lbound;

%  plot(params(:,1), 'b')

  [dm,ds] = mymean(diff(params, [], 1), 1);

  imf = emdc([], params(:,1), true, 3);
%  p0(:,1) = sum(imf(end-3:end,:), 1).';
  p0(:,1) = imf(end,:).';
  bound(1) = ds(1)/2;

  imf = emdc([], params(:,2), true, 4);
  %p0(:,2) = sum(imf(end-5:end,:), 1).';
  p0(:,2) = imf(end,:).';
  bound(2) = ds(2)/4;

  [m,s] = mymean(params(:,3));
  p0(:,3) = m;
  bound(3) = s/10;

  lbound = bsxfun(@minus, p0, bound);
  ubound = bsxfun(@plus, p0, bound);

  p0(p0(:,1)<0,1)=0;
  lbound(lbound(:,1)<0,1)=0;

  [niter, params] = fit_cminpack_mex(p0.', values.', coef, lbound.', ubound.', dpos, smoothed.');
  
  gaussians = params.';
  goods = (gaussians(:,3) ~= -1);
  gaussians(~goods, :) = NaN;
  gaussians(gaussians(:,1) < 0, 1) = 0;
%  plot(gaussians(:,1), 'c')

%  keyboard

  return;
end
