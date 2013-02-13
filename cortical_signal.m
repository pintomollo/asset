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
  for i=1:ninit
    nimg = randi(nframes);
    ph = double(load_data(mymovie.cortex, nimg));
    cortex = mymovie.data.cortex(nimg).carth;

    if (~isempty(cortex))
      [priors(i,:)] = estimate_peak(ph, cortex, opts);
    end
  end
  
  priors = median(priors, 1);
%  keyboard

  for i=1:nframes
    %bounds = [-Inf bin_step 0 0 -Inf 0 0; ...
    %           Inf Inf Inf Inf Inf Inf Inf];

    %bounds_invag = [-Inf bin_step -Inf -Inf 0 0; ...
    %                 Inf Inf Inf Inf Inf Inf];

    nimg = i;
    %nimg = randi(nframes)
    %nimg = i + 213
    %nimg = 80

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
      figure;imagesc(img);

      safety = ceil(2 / opts.pixel_size);
      path = draw_ellipse(mymovie.data.centers(:, nimg), mymovie.data.axes_length(:, nimg)+safety, mymovie.data.orientations(nimg));
      mask = roipoly(img, path(:,1), path(:,2));
      img2 = img;
      img2(mask) = NaN;
      vertical_bkg = nanmean(img2,2);
      img = bsxfun(@minus, img, vertical_bkg);
      figure;imagesc(img);
      img2 = img;
      img2(mask) = NaN;
      horizontal_bkg = nanmean(img2,1);
      img = bsxfun(@minus, img, horizontal_bkg);
      img(img < 0) = 0;

      figure;imagesc(img);
      keyboard

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
      [gs, dperp] = estimate_signal(img, cortex, priors, opts);

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

  if (opts.recompute|~isfield(mymovie.data, 'domain')|isempty(mymovie.data.domain))
    mymovie = carth2normalized(mymovie, opts);

    params = opts.quantification.params;
    weights = opts.quantification.weights;
    init_params = opts.quantification.init_params;

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

function [range, center] = get_peak(x, y)

  % Get only the peak
  dy = differentiator(x, y);

  ymax = find(dy(1:end-1) > 0 & dy(2:end) <= 0)+1;
  if (isempty(ymax))
    [junk, center] = min(abs(x));
    center = x(center);
  else
    [junk, center] = min(abs(x(ymax)));
    center = x(ymax(center));
  end

  x = x - center;

  mins = find(dy(1:end-1) < 0 & dy(2:end) >= 0)+1;
  maxs = [];
  lmin = mins(x(mins) < 0);
  if (isempty(lmin))
    maxs = find(dy(1:end-2) > dy(2:end-1) & dy(2:end-1) <= dy(3:end)) + 1;
    lmin = maxs(x(maxs) < 0);
    if (isempty(lmin))
      lmin = 1;
    end
  end
  rmin = mins(x(mins) > 0);
  if (isempty(rmin))
    if (isempty(maxs))
      maxs = find(dy(1:end-2) > dy(2:end-1) & dy(2:end-1) <= dy(3:end)) + 1;
    end
    rmin = maxs(x(maxs) > 0);
    if (isempty(rmin))
      rmin = length(x);
    end
  end

  corrs = ones(size(lmin));
  all_indx = zeros(size(lmin));
  for i=1:length(lmin)
    [junk, all_indx(i)] = min(abs(x(rmin) + x(lmin(i))));
    indexes = [lmin(i):rmin(all_indx(i))];
    corrs(i) = gaussian_correlation(x(indexes), y(indexes), 0.05);
  end

  [junk, best_min] = max(corrs);

  if (isempty(best_min))
    best_min = 1;
  end

  range = mean([-x(lmin(best_min)), x(rmin(all_indx(best_min)))]);

  % Approx 1/20
  center = [center range/2.45];

  range = (x >= -range & x <= range);

  return;
end

function [sigma, ks] = normality_test(x, y, factor)
% Kolmogorov-Smirnov test

  if (nargin == 2)
    factor = 3;
  end

  min_length = 2;
  half = find(x==0);
  max_length = floor(min(half-1, length(x)-half)/(2*factor)) - 1;
  ks = NaN(max_length, 1);

  if (max_length < min_length)
    sigma = NaN;
    ks = NaN;
    return;
  end

  for i=min_length:max_length
    w = x(half+i);
    pos = (abs(x) <= factor*w);
    bkg = interp1(x(~pos), y(~pos), x(pos), 'linear');

    tmp_x = x(pos) / w;
    tmp_y = y(pos) - bkg;

    tmp_y = cumsum(tmp_y);
    tmp_y = tmp_y / tmp_y(end);

    norm_vals = cdf('norm', tmp_x, 0, 1);
    delta = abs(norm_vals - tmp_y);
    ks(i, 1) = max(delta);
  end

  [ks, indx] = min(ks, [], 1);
  sigma = x(half+indx);

  return;
end

function correl = gaussian_correlation(x, y, thresh)

  y = y - mean(y([1 end]));
  ampl = max(y) / (1 - thresh);
  sigma = mean(real(sqrt(-(x([1 end]).^2)/(2*log(thresh)))));

  vals = ampl*(exp(-(x.^2) / (2*sigma^2)) - thresh);
  correl = corrcoef(y(:), vals(:));
  correl = correl(1,2);

  return;
end

function [peak_params, dperp, dpos] = estimate_peak(ph, cortex, priors, opts)

  if (nargin == 3)
    opts = priors;
    priors = [];
  end

  ph = gaussian_mex(ph, 1.5);
  [ph_values, dperp, dpos] = perpendicular_sampling(ph, cortex, opts);

  y = mymean(ph_values, 1);
  x = dpos;

  valids = isfinite(y);

  y = y(valids);
  x = x(valids);

%  imf = emdc(x, y);
%  if (size(imf, 1) > 2)
%    smoothed = sum(imf(end-1:end, :), 1);
%  else
%    smoothed = imf(end, :);
%  end

  dy = differentiator(x, y, true);
  ddy = -differentiator(x, dy, true);
  dddy = differentiator(x, ddy, true);

  maxs = find(dddy(1:end-1) > 0 & dddy(2:end) <= 0);
  ks = NaN(length(maxs), 2);

  for i=1:length(maxs)
    [ks(i,1), ks(i,2)] = normality_test(x - x(maxs(i)), y);
  end

  if (~isempty(priors))
    ks(ks(:,1) < priors(2)-1 | ks(:,1) > priors(2)+1, :) = NaN;
  end

  [junk, indx] = min(ks(:,2));
  peak_params = [x(maxs(indx)) ks(indx, 1)];

  if (~isempty(priors) & ~isfinite(peak_params(1)))
    peak_params = priors;
  end

  %keyboard

  return;

  [junk, indx] = min(abs(x(maxs)));
  indx = maxs(indx);

  if (isempty(indx))
    center = 0;
  else
    center = x(indx);
  end

  [peak, peak_params] = get_peak(x - center, ddy);

  peak_params(1) = center + peak_params(1);
  peak_params(2) = sqrt(3)*peak_params(2);

  
%  ampls(ampls < 0) = 0;
%  p0 = [ampls(:) repmat(peak_prior, npts, 1)];

  lbound = [0 -0.25*peak_params(2)+peak_params(1) peak_params(2)];
  ubound = [max(smoothed(:)) 0.25*peak_params(2)+peak_params(1) peak_params(2)*3];

  test = y;
  test(peak) = interp1(x(~peak), y(~peak), x(peak), 'linear');

  keyboard

  %[niter, params] = fit_gaussian_mex(p0.', values.', 400, 1e-6, lbound.', ubound.', dpos, smoothed.');
  [niter, params] = fit_cminpack_mex(p0.', values.', 400, 1e-6, lbound.', ubound.', dpos, smoothed.');



  return;
end

%function [gaussians, dperp] = estimate_signal(img, cortex, dperp, dpos, peak_prior, opts)
function [gaussians, dperp] = estimate_signal(img, cortex, peak_prior, opts)
 
  [values, dperp, dpos] = perpendicular_sampling(img, cortex, opts);
  
  halfed = gaussian_mex(img, peak_prior(2)/2);
  [halfed] = perpendicular_sampling(halfed, cortex, dperp, dpos, opts);

  %{
  img = gaussian_mex(img, peak_prior(2));
  [smoothed] = perpendicular_sampling(img, cortex, dperp, dpos, opts);

  g = exp(-(dpos-peak_prior(1)).^2 / (2*peak_prior(2)^2));

  [X,Y] = meshgrid(dpos, [1:size(values,1)].');

  bkg = smoothed;
  coef = 3;

  % Approx 1/5
  peak_dist = coef*peak_prior(2);

  gauss = (abs(dpos-peak_prior(1)) < peak_dist);

  prc = [0 0.25 0.5 0.75 1];
  vals = prctile(halfed(:,gauss).', prc*100);

  coefs = [prc.' ones(length(prc), 1)] \ vals;
  ampls = coefs(1,:).';

  imf = emdc([], ampls, true);
  if (size(imf, 1) > 4)
    ampls = sum(imf(end-3:end, :), 1);
  else
    ampls = imf(end, :);
  end

  ampls(ampls < 0) = 0;
  bkg = halfed;
  bkg(:,gauss) = interp2(X(:,~gauss), Y(:,~gauss), bkg(:,~gauss), X(:,gauss), Y(:, gauss), 'linear');
  approx1 = bkg + bsxfun(@times, ampls.', g);

  gaussians = perform_fit(ampls, peak_prior, dpos, values, halfed, coef);

  p = bsxfun(@minus, dpos, gaussians(:,2));
  gauss = logical(median(double(bsxfun(@lt, abs(p), gaussians(:,3)*coef)), 1));

  bkg = halfed;
  bkg(:,gauss) = interp2(X(:,~gauss), Y(:,~gauss), bkg(:,~gauss), X(:,gauss), Y(:,gauss), 'linear');
  final1 = bkg + bsxfun(@times, gaussians(:,1), exp(bsxfun(@rdivide, -p.^2, 2*gaussians(:,3).^2)));

  %------------------------------------------

  coef = 3;
  peak_dist = coef*peak_prior(2);

  gauss2 = (abs(dpos-peak_prior(1)) < peak_dist);

  prc = [0 0.25 0.5 0.75 1];
  vals = prctile(smoothed(:,gauss2).', prc*100);

  coefs = [prc.' ones(length(prc), 1)] \ vals;
  ampls = coefs(1,:).';

  imf = emdc([], ampls, true);
  if (size(imf, 1) > 4)
    ampls = sum(imf(end-3:end, :), 1);
  else
    ampls = imf(end, :);
  end

  ampls(ampls < 0) = 0;
  bkg = smoothed;
  bkg(:,gauss2) = interp2(X(:,~gauss2), Y(:,~gauss2), bkg(:,~gauss2), X(:,gauss2), Y(:, gauss2), 'linear');
  approx2 = bkg + bsxfun(@times, ampls.', g);

  gaussians = perform_fit(ampls, peak_prior, dpos, values, smoothed, coef);

  p = bsxfun(@minus, dpos, gaussians(:,2));
  gauss = logical(median(double(bsxfun(@lt, abs(p), gaussians(:,3)*coef)), 1));

  bkg = smoothed;
  bkg(:,gauss) = interp2(X(:,~gauss), Y(:,~gauss), bkg(:,~gauss), X(:,gauss), Y(:,gauss), 'linear');
  final2 = bkg + bsxfun(@times, gaussians(:,1), exp(bsxfun(@rdivide, -p.^2, 2*gaussians(:,3).^2)));
  %}

  %------------------------------------------

  coef = 3;
  peak_dist = coef*peak_prior(2);

  gauss3 = (abs(dpos-peak_prior(1)) < peak_dist);

  [X,Y] = meshgrid(dpos, [1:size(values,1)].');

  bkg = halfed;
  bkg(:,gauss3) = interp2(X(:,~gauss3), Y(:,~gauss3), bkg(:,~gauss3), X(:,gauss3), Y(:, gauss3), 'linear');

  s = halfed - bkg;
  [maxs, indxs] = max(s(:, gauss3), [], 2);
  avgs = mymean(s(:, gauss3), 2);
  % Average over [-3 3]*sigma of a gaussian is ~0.42*ampl
  ampls = (maxs-avgs)/0.58;

  imf = emdc([], ampls, true);
  if (size(imf, 1) > 4)
    ampls = sum(imf(end-3:end, :), 1);
  else
    ampls = imf(end, :);
  end

  %bkg = halfed;
  %bkg(:,gauss3) = interp2(X(:,~gauss3), Y(:,~gauss3), bkg(:,~gauss3), X(:,gauss3), Y(:, gauss3), 'linear');

  ampls(ampls < 0) = 0;
  %approx3 = bkg + bsxfun(@times, ampls.', g);

  gaussians = perform_fit(ampls, peak_prior, dpos, values, halfed, coef);

  %p = bsxfun(@minus, dpos, gaussians(:,2));
  %gauss = logical(median(double(bsxfun(@lt, abs(p), gaussians(:,3)*coef)), 1));

  %bkg = halfed;
  %bkg(:,gauss) = interp2(X(:,~gauss), Y(:,~gauss), bkg(:,~gauss), X(:,gauss), Y(:,gauss), 'linear');
  %final3 = bkg + bsxfun(@times, gaussians(:,1), exp(bsxfun(@rdivide, -p.^2, 2*gaussians(:,3).^2)));

  %------------------------------------------
  %{
  coef = 3;
  peak_dist = coef*peak_prior(2);

  gauss4 = (abs(dpos-peak_prior(1)) < peak_dist);

  bkg = smoothed;
  bkg(:,gauss4) = interp2(X(:,~gauss4), Y(:,~gauss4), bkg(:,~gauss4), X(:,gauss4), Y(:, gauss4), 'linear');

  s = smoothed - bkg;
  [maxs, indxs] = max(s(:, gauss4), [], 2);
  avgs = mymean(s(:, gauss4), 2);
  % Average over [-3 3]*sigma of a gaussian is ~0.42*ampl
  ampls = (maxs-avgs)/0.58;

  imf = emdc([], ampls, true);
  if (size(imf, 1) > 4)
    ampls = sum(imf(end-3:end, :), 1);
  else
    ampls = imf(end, :);
  end

  bkg = smoothed;
  bkg(:,gauss4) = interp2(X(:,~gauss4), Y(:,~gauss4), bkg(:,~gauss4), X(:,gauss4), Y(:, gauss4), 'linear');

  ampls(ampls < 0) = 0;
  approx4 = bkg + bsxfun(@times, ampls.', g);

  gaussians = perform_fit(ampls, peak_prior, dpos, values, smoothed, coef);

  p = bsxfun(@minus, dpos, gaussians(:,2));
  gauss = logical(median(double(bsxfun(@lt, abs(p), gaussians(:,3)*coef)), 1));

  bkg = smoothed;
  bkg(:,gauss) = interp2(X(:,~gauss), Y(:,~gauss), bkg(:,~gauss), X(:,gauss), Y(:,gauss), 'linear');
  final4 = bkg + bsxfun(@times, gaussians(:,1), exp(bsxfun(@rdivide, -p.^2, 2*gaussians(:,3).^2)));

  [[nanmean(nanmean(abs(values - approx1))); ...
  nanmean(nanmean(abs(values - approx2))); ...
  nanmean(nanmean(abs(values - approx3))); ...
  nanmean(nanmean(abs(values - approx4)))], ...
  [nanmean(nanmean(abs(values - final1))); ...
  nanmean(nanmean(abs(values - final2))); ...
  nanmean(nanmean(abs(values - final3))); ...
  nanmean(nanmean(abs(values - final4)))]]
  %}
  
  return;
end

function gaussians = perform_fit(ampls, peak_prior, dpos, values, smoothed, coef)
%  [maxs, indxs] = max(smoothed(:, gauss), [], 2);
%  avgs = mymean(smoothed(:, gauss), 2);
  % Average over [-1.79 1.79]*sigma of a gaussian is ~0.65*ampl
%  ampls = (maxs-avgs)/0.35;

  npts = length(ampls);

  p0 = [ampls(:) repmat(peak_prior, npts, 1)];
  %lbound = [0 -2*peak_prior(2)+peak_prior(1) peak_prior(2)*0.75];
  %ubound = [max(smoothed(:)) 2*peak_prior(2)+peak_prior(1) 2*peak_prior(2)];

  lbound = [0 -2*peak_prior(2)+peak_prior(1) peak_prior(2)*0.95];
  ubound = [max(smoothed(:)) 2*peak_prior(2)+peak_prior(1) peak_prior(2)*1.05];

  %[niter, params] = fit_gaussian_mex(p0.', values.', 400, 1e-6, lbound.', ubound.', dpos, smoothed.');
  [niter, params] = fit_cminpack_mex(p0.', values.', coef, 1e-6, lbound.', ubound.', dpos, smoothed.');

  bound = lbound;
  params = params.';

%  plot(params(:,1), 'b')

  [dm,ds] = mymean(diff(params, [], 1), 1);

  half = round(npts/2);

  imf = emdc([], params(:,1), true);
  p0(:,1) = sum(imf(end-3:end,:), 1).';
  bound(1) = ds(1)/2;

  imf = emdc([], params(:,2), true);
  p0(:,2) = sum(imf(end-2:end,:), 1).';
  bound(2) = ds(2)/4;

  [m,s] = mymean(params(:,3));
  p0(:,3) = m;
  bound(3) = s/10;

  lbound = bsxfun(@minus, p0, bound);
  ubound = bsxfun(@plus, p0, bound);

  p0(p0(:,1)<0,1)=0;
  lbound(lbound(:,1)<0,1)=0;

  %[niter, params] = fit_gaussian_mex(p0.', values.', 400, 1e-6, lbound.', ubound.', dpos, smoothed.');
  [niter, params] = fit_cminpack_mex(p0.', values.', coef, 1e-6, lbound.', ubound.', dpos, smoothed.');
  
  gaussians = params.';
%  plot(gaussians(:,1), 'c')

%  keyboard

  return;
end
