function mymovie = cortical_signal(mymovie, opts)

  type = opts.segmentation_type;
  if (strncmp(type, 'markers', 7))
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

  if (opts.recompute | ~isfield(mymovie.data, 'eggshell') | isempty(mymovie.data.eggshell))
    mymovie = duplicate_segmentation(mymovie, 'data', opts);
  end

  if (opts.quantification.use_ruffles)
    if (~isfield(mymovie.(type), 'ruffles') | isempty(mymovie.(type).ruffles))
      mymovie = find_ruffles(mymovie, opts);
      mymovie = follow_invaginations(mymovie, opts);
    elseif (~isfield(mymovie.(type).ruffles, 'paths') | isempty(mymovie.(type).ruffles(1).paths))
      do_it = true;
      for j=1:nframes
        if (~isempty(mymovie.(type).ruffles(j).paths))
          do_it = false;
        end
      end
      if (do_it)
        mymovie = follow_invaginations(mymovie, opts);
      end
    end
  end

  bounds = [-Inf 0 0 0 -Inf 0 0; ...
             Inf Inf Inf Inf Inf Inf Inf];

  bounds_invag = [-Inf 0 -Inf -Inf 0 0; ...
                   Inf Inf Inf Inf Inf Inf];

  optims = optimset('Display', 'off');
  
  for i=1:nframes
    %nimg = randi(nframes)
    %nimg = i;
    %nimg = i + 85
    nimg = 91
    %keyboard

    %keyboard

    cortex = mymovie.data.cortex(nimg).carth;

    if (opts.recompute||(~isfield(mymovie.data, 'quantification'))||(length(mymovie.data.quantification) < nimg)||~isfield(mymovie.data.quantification(nimg), 'front')||isempty(mymovie.data.quantification(nimg).front))
      img = double(load_data(mymovie.data, nimg));
      img = mask_neighbors(img, mymovie.data.centers(:,nimg), mymovie.data.axes_length(:,nimg), mymovie.data.orientations(1,nimg), mymovie.data.neighbors(nimg), opts);

      ell_cortex = carth2elliptic(cortex, mymovie.data.centers(:, nimg), mymovie.data.axes_length(:, nimg), mymovie.data.orientations(1, nimg));

      ell_pos = ell_cortex(:,1);
      ell_pos = unique([ell_pos; 0; pi]);
      ell_pos = ell_pos([diff([ell_pos; ell_pos(1)+2*pi]) > 1e-10]);

      ell_cortex = interp_elliptic(ell_cortex, ell_pos);

      cortex = elliptic2carth(ell_cortex, mymovie.data.centers(:, nimg), mymovie.data.axes_length(:, nimg), mymovie.data.orientations(1, nimg));

      if (opts.quantification.use_ruffles)
        [ cortex, rescale] = insert_ruffles(cortex, mymovie.(type).ruffles(nimg).paths);
      else
        rescale = false(size(cortex, 1), 1);
      end

      tmp_cortex = [cortex(end,:); cortex; cortex(1,:)];
      dpos = (tmp_cortex(3:end, :) - tmp_cortex(1:end-2, :)) / 2;
      dpts = diff(cortex(1:end-1, :));

      dperp = [-dpos(:,2) dpos(:,1)];
      nulls = all(dperp == 0, 2);
      dperp(nulls, :) = dpts(nulls, :);

      dperp = bsxfun(@rdivide, dperp, hypot(dperp(:,1), dperp(:,2))); 
      nbins = 64;

      dpos = [-nbins:nbins];
      nbins = length(dpos);

      all_pos_x = bsxfun(@plus, dperp(:,1) * dpos, cortex(:,1));
      all_pos_y = bsxfun(@plus, dperp(:,2) * dpos, cortex(:,2));

      all_pos_x = all_pos_x(:);
      all_pos_y = all_pos_y(:);

      values = bilinear(img, all_pos_x, all_pos_y);
      values = reshape(values, [size(cortex,1) nbins]);

      npts = size(values, 1);
      slopes = NaN(npts, 2);

      [params, bounds, slopes(~rescale, 1)] = estimate_mean(dpos, values(~rescale, :), false, bounds, optims);
      bounds_invag(:, 1) = bounds(:, 1) + bounds(:, 2);
      bounds_invag(:, 2) = 2*bounds(:,2);
      [params_invag, bounds_invag, slopes(rescale, :)] = estimate_mean(dpos, values(rescale, :), true, bounds_invag, optims);

      %projection = mymean(values(~rescale, :), 1);
      %valids = ~isnan(projection);

      %valid_projection = projection(valids);
      %valid_dpos = dpos(valids);

      %imf = emdc(valid_dpos, valid_projection);
      %valid_projection = sum(imf(2:end, :), 1);

      %keyboard

      %params = estimate_front(valid_dpos, valid_projection, []);
      %best1 = nlinfit(valid_dpos, valid_projection, @front, params);
      %best2 = lsqcurvefit(@front, params, valid_dpos, valid_projection, lbound, ubound, optims);
      
      %obj_function = valid_projection;
      %position = valid_dpos;
      %best2 = lsqnonlin(@fit_front, params);

      %keyboard
      
      %figure;plot(valid_dpos, valid_projection, 'b');
      %hold on;
      %projection = mymean(values(rescale, :), 1);
      %valids = ~isnan(projection);
      %plot(dpos(valids), projection(valids), 'r');

      %figure;imshow(imnorm(img));
      %hold on;
      %myplot(cortex)
      %myplot(cortex(rescale, :), 'r');

      %figure;imagesc(values(rescale,:));

      %keyboard

      %continue;
      %plot(dpos, front(params, dpos), 'c');
      %plot(dpos, front(best1, dpos), 'g');
      %plot(dpos, front(best2, dpos), 'r');

      %projection = -mymean(diff(values(~rescale, valid_dpos < best2(1) - 3*best2(2)), [], 2), 2);
      %imf = emdc([1:length(projection)], projection);
      %slopes = NaN(size(values, 1), 1);
      %slopes(~rescale) = imf(end, :);

      %figure;subplot(121);plot(projection);hold on;plot(imf(end, :), 'r');
      %subplot(122);imagesc(values);
     %keyboard
      %continue
      %params = estimate_front(valid_dpos, valid_projection, params);
      %best2 = nlinfit(valid_dpos, projection(valids), @front, params);


      %params1 = estimate_sigmoid(valid_dpos, valid_projection);
      %estim1 = sigmoid(params1, valid_dpos);

      %estim5 = piecewise_background(best2, valid_dpos);
      %peak_range = get_peak(valid_dpos, valid_projection-estim5);

      %params2 = estimate_gaussian(valid_dpos(peak_range), valid_projection(peak_range)-estim1(peak_range));
      %estim2 = gaussian([params2(1:3) 0], valid_dpos);

      %params3 = estimate_sigmoid(valid_dpos, valid_projection - estim2);
      %estim3 = sigmoid(params3, valid_dpos);

      %peak_range = get_peak(valid_dpos, valid_projection - estim3);
      %params4 = estimate_gaussian(valid_dpos(peak_range), valid_projection(peak_range)-estim3(peak_range));

      %params5 = params3;
      %params5(4) = params5(4) + params4(4);
      %estim5 = sigmoid(params5, dpos);

      %peak_range = correct_peak(peak_range, valids);

      signal = NaN(npts, 7); 
      %signal2 = NaN(npts, 7); 
      %smoothed = NaN(size(values));

      for j=1:npts

        if (rescale(j))
        %      continue;
          curr_bounds = bounds_invag;
          curr_params = [params_invag(1:2) slopes(j, :) params_invag(end-1:end)];
        else
          curr_bounds = bounds;
          curr_params = [params(1:4) slopes(j, 1) params(end-1:end)];
        end
          
        line = values(j,:);
        valids = ~isnan(line);

        if (~any(valids))
          %all_pos_x = reshape(all_pos_x, [size(cortex,1) nbins]);
          %all_pos_y = reshape(all_pos_y, [size(cortex,1) nbins]);

          %[all_pos_x(j, :); all_pos_y(j, :)]

          %j

          warning(['No valid data point for quantification in frame ' num2str(nimg) ', line ' num2str(j) '.'])

          continue;
        end

        valid_line = line(valids);
        valid_dpos = dpos(valids);
        %valid_range = valids & peak_range;

        tmp_val = emdc(valid_dpos, valid_line);
        smoothed = sum(tmp_val(end-1:end, :), 1);
        
        %p1 = estimate_front(dpos(valids), line(valids), curr_params);
        p2 = estimate_front(valid_dpos, smoothed, curr_params, rescale(j));
        
        %estim = piecewise_background(curr_params, dpos(valid_range));
        %signal(j, :) = estimate_gaussian(dpos(valid_range), smoothed(j, valid_range)-estim);
        %tmp_params = estimate_gaussian(dpos(valid_range), smoothed(j, valid_range)-estim);
        %tmp_params = [tmp_params(1:2) curr_params(3:5) curr_params(6)+tmp_params(4) tmp_params(3)];

        try
        %signal(j, :) = nlinfit(dpos(valids), valid_line, @front, tmp_params);
          %s1 = nlinfit(dpos(valids).', valid_line.', @front, p1, statset('Robust', 'on'));
          s1 = lsqcurvefit(@front, p2, valid_dpos.', valid_line.', curr_bounds(1,:), curr_bounds(2, :), optims);
          %s2 = nlinfit(dpos(valids).', valid_line.', @front, p2, statset('Robust', 'on'));
        catch
          keyboard
        end

        if (any(front(s1, dpos) < 0))
          keyboard
        %elseif (any(front(s2, dpos) < 0))
        %  keyboard
        end

        if (rescale(j))
          s1 = [s1(1:2) NaN s1(3:end)];
        end

        signal(j, :) = s1;
        %signal2(j, :) = s2;

        %figure;plot(dpos(valids), valid_line, 'b');
        %hold on;
        %plot(dpos, front(p1, dpos), 'm');
        %plot(dpos, front(p2, dpos), 'c');
        %plot(dpos, front(s1, dpos), 'r');
        %plot(dpos, front(s2, dpos), 'g');

        %keyboard

        %if (rescale(j))
        %  signal(j, 3) = signal(j, 3) / 2;
        %end
      end

      mymovie.data.quantification(nimg).front = signal(:, [1 2 end]);
      mymovie.data.quantification(nimg).bkg = signal(:, [3:end-1]);
      mymovie.data.quantification(nimg).carth = cortex;

      %figure;imagesc(gaussian(signal(:, [1 2 7]), dpos));

      %if (false)
      %figure;imshow(imnorm(img));hold on;
      %quiver(cortex(:,1), cortex(:, 2), dperp(:, 1), dperp(:, 2), 'r')

      %estim4 = gaussian([params4(1:3) 0], dpos);
      %clim = [min(values(:)) max(values(:))];
      %figure;imagesc(gaussian(signal, dpos));
      %figure;imagesc(bsxfun(@plus, gaussian(signal, dpos), estim5), clim);
      %figure;imagesc(values, clim);
      %figure;imagesc(smoothed);
      %figure;plot(dpos, projection);hold on;plot(dpos(peak_range), projection(peak_range), 'r');
      %plot(dpos, estim5, 'k');plot(dpos, estim4, 'g');
      %plot(dpos, estim5 + estim4, 'm')
      %plot(dpos, estim3, 'm');plot(dpos, estim2, 'y');
      %plot(dpos, estim1, 'c');
      %keyboard
      %end
      %close all
    else

      tmp_cortex = [cortex(end,:); cortex; cortex(1,:)];
      dpos = (tmp_cortex(3:end, :) - tmp_cortex(1:end-2, :)) / 2;

      dperp = [-dpos(:,2) dpos(:,1)];
      dperp = bsxfun(@rdivide, dperp, hypot(dperp(:,1), dperp(:,2))); 

      signal = mymovie.data.quantification(nimg).front;
    end
    
    if (opts.recompute|length(mymovie.data.quantification) < nimg|~isfield(mymovie.data.quantification(nimg), 'cortex')|isempty(mymovie.data.quantification(nimg).cortex))
      values = extract_ridge(signal, cortex, dperp, opts);
      mymovie.data.quantification(nimg).cortex = values;
    end
  end

  %if (false)
  if (opts.recompute|~isfield(mymovie.data, 'domain')|isempty(mymovie.data.domain))
    mymovie = carth2normalized(mymovie, opts);

  params = opts.quantification.params;
  weights = opts.quantification.weights;
  init_params = opts.quantification.init_params;

    img = gather_quantification(mymovie, opts);
    img = imnorm(img);
    %opts.quantification.params.init = sub2ind(size(img,2)*[1 1], [1:size(img,2)],[1:size(img,2)]);
      path = dynamic_prog_2d(img, params, @weight_domain, weights, @init_domain, init_params, opts);

    %path = dynamic_prog_2d(img, opts.quantification.params, @weight_domain_borders, opts.quantification.weights, opts);
    mymovie.data.domain = path / size(img, 2);
  end

  return;

%  function err = fit_front(current_parameters)
%
%    lparams = (current_parameters < lbound);
%    uparams = (current_parameters > ubound);
%
%    val = front(current_parameters, position);
%    err = abs(val - obj_function);
%    err = mean(err) + std(err) + exp(10*(sum(abs(current_parameters(lparams) - lbound(lparams))) + sum(abs(current_parameters(uparams) - ubound(uparams)))));
%
%    return;
%  end
end

function [params, bounds, slopes] = estimate_mean(x, ys, is_invag, bounds, optims)

  if (isempty(ys))  
    params = [];
    bounds = [];
    slopes = zeros(0,2);

    return;
  end

  y = mymean(ys, 1);
  valids = isfinite(y);

  y = y(valids);
  x = x(valids);

  imf = emdc(x, y);
  y = sum(imf(2:end, :), 1);

  params = estimate_front(x, y, [], is_invag);
  params = lsqcurvefit(@front, params, x, y, bounds(1,:), bounds(2, :), optims);

  y = -mymean(diff(ys(:, x < params(1) - 3*params(2)), [], 2), 2);

  if (isempty(y))
    slopes = ones(size(ys, 1), 1) * params(end-1);
  else
    imf = emdc([], y);
    slopes = imf(end, :).';
  end

  if (is_invag)
    y = -mymean(diff(ys(:, x > params(1) + 3*params(2)), [], 2), 2);

    if (isempty(y))
      slopes = [slopes ones(size(slopes)) * params(end-2)];
    else
      imf = emdc([], y);
      slopes = [slopes, imf(end, :).'];
    end
  end

  bounds(1, 1) = params(1) - 2*params(2);
  bounds(2, 1) = params(1) + 2*params(2);
  bounds(2, 2) = 2.5*params(2);

  return;
end

function new_peak = correct_peak(peak, valids)

  new_peak = false(size(valids));
  new_peak(valids) = peak;

  start_indx = find(new_peak, 1, 'first');
  end_indx = find(new_peak, 1, 'last');

  new_peak(start_indx:end_indx) = true;

  return;
end

function new_params = estimate_front(x, y, params, is_invag)

  params_sigmoid = estimate_piecewise(x, y, params, is_invag);

  %estim = sigmoid(params_sigmoid, x);
  estim = piecewise_background(params_sigmoid, x);
  peak_range = get_peak(x, y-estim);
  params_gaussian = estimate_gaussian(x(peak_range), y(peak_range)-estim(peak_range));

  %new_params = [params_sigmoid(1:end-1) params_sigmoid(end)+params_gaussian(end) params_gaussian(1:end-1)];
  new_params = [params_sigmoid(1) params_gaussian(2) params_sigmoid(3:end-1) params_sigmoid(end)+params_gaussian(end) params_gaussian(3:end-1)];

  return;
end

function params = estimate_sigmoid(x, y)

  % Derivative of sigmoid is a gaussian, which we know how to estimate !
  dy = -diff(y);
  pos = x(2:end);

  p = estimate_gaussian(pos, dy);
  estim = sigmoid([p(1:2) 1 0], x);

  vals = [estim.' ones(size(x)).'] \ y.';
  params = [p(1:2) vals.'];

  return;
end

function [range, center] = get_peak(x, y)

  %keyboard

  % Get only the peak
  dy = differentiator(x, y);

  ymax = find(dy(1:end-1) > 0 & dy(2:end) <= 0)+1;
  if (isempty(ymax))
    [~, center] = min(abs(x));
    center = x(center);
  else
    [~, center] = min(abs(x(ymax)));
    center = x(ymax(center));
  end

  %figure;plot(x,y);
  %hold on;
  %plot(x,dy, 'k')
  %keyboard

  x = x - center;

  %lmin = find(x(2:end-1) < 0 & dy(1:end-1) < 0 & dy(2:end) >= 0)+1;
  %rmin = find(x(2:end-1) > 0 & dy(1:end-1) <= 0 & dy(2:end) > 0)+1;

  %mins = find(dy(1:end-2) < dy(2:end-1) & dy(2:end-1) >= dy(3:end)) + 2;
  mins = find(y(1:end-2) > y(2:end-1) & y(2:end-1) <= y(3:end)) + 1;
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
    [~, all_indx(i)] = min(abs(x(rmin) + x(lmin(i))));
    indexes = [lmin(i):rmin(all_indx(i))];
    corrs(i) = gaussian_correlation(x(indexes), y(indexes), 0.05);
  end

  [~, best_min] = max(corrs);

  %width = abs(bsxfun(@minus, x(rmin), x(lmin).'));
  %dist = abs(bsxfun(@plus, x(rmin), x(lmin).'));
  %dist(width < 30) = Inf;
  %[dist, indexes] = min(dist, [], 2);
  %best_min = find(dist(2:end) <= dist(1:end-1), 1, 'last') + 1;
  if (isempty(best_min))
    best_min = 1;
    %best_min = length(dist);
  end

  range = min(-x(lmin(best_min)), x(rmin(all_indx(best_min))));
  center = [center range];
  range = (x >= -range & x <= range);

  return;
end

function correl = gaussian_correlation(x, y, thresh)

  y = y - mean(y([1 end]));
  ampl = max(y) / (1 - thresh);
  sigma = mean(real(sqrt(-(x([1 end]).^2)/(2*log(thresh)))));

  vals = ampl*(exp(-(x.^2) / (2*sigma^2)) - thresh);
  correl = corrcoef(y(:), vals(:));
  correl = correl(1,2);

  %figure;plot(x,y,'b');
  %hold on;
  %plot(x,vals,'r')

  return;
end

function params = estimate_gaussian(x,y)

  orig_y = y;
  y = y - min(y);
  
  y = y / sum(y);
  center = sum(x .* y);

  x = (x - center).^2;
  sigma = sum(y .* x);

  vals = [exp(-x / (2*sigma)).' ones(size(x)).'] \ orig_y.';
  [ampl, bkg] = deal(vals(1), vals(2));
  if (ampl < 0)
    ampl = 0;
  end
  sigma = sqrt(sigma);

  params = [center sigma ampl bkg];

  return;
end

function y = front(p, x)

  fliped = false;
  if (size(x,1) ~= 1)
    x = x.';
    fliped = true;
  end

  if (numel(p) == 6)
    p = [p(1:2) NaN p(3:end)];
  end

  %y = sigmoid(p(1:4), x) + gaussian([p(5:end) 0], x);
  %y = piecewise_background(p(1:6), x) + gaussian([p(end-2:end) 0], x);
  y = piecewise_background(p(:, 1:end-1), x) + gaussian([p(:, [1:2 end]) zeros(size(p(:,1)))], x);

  if (fliped)
    y = y.';
  end

  return;
end

function y = sigmoid(p, x)

  y = bsxfun(@plus, bsxfun(@times, ((1 + erf(bsxfun(@rdivide, bsxfun(@plus, -x, p(:, 1)), p(:, 2)*sqrt(2)))) ./ 2), p(:, 3)), p(:, 4));

  return;
end

function y = gaussian(p, x)

  if (any(isnan(p)) | any(isinf(p)))
    y = zeros(size(x));

    return;

  elseif (size(p, 2) < 4)
    p = [p zeros(size(p, 1), 1)];
  end

  y = bsxfun(@plus, bsxfun(@times, p(:, 3), exp(bsxfun(@rdivide, -(bsxfun(@minus,x,p(:,1)).^2), (2*(p(:, 2).^2))))), p(:, 4));

  return;
end

function params = estimate_piecewise(x, y, params, is_invag)

  if (nargin < 3)
    params = [];
    is_invag = false;
  elseif (nargin == 3)
    if (islogical(params))
      is_invag = params;
      params = [];
    else
      is_invag = false;
    end
  end

  dy = differentiator(x, y, true);

  if (~isempty(params))
    gauss = imnorm(gaussian([params(1) 2*params(2) params(end) 0.1 * params(end)], x));
    ddy = -differentiator(x, dy .* gauss, true);
  else
    ddy = -differentiator(x, dy, true);
  end

  %dy = diff(y) ./ diff(x);
  %ddy = -diff(dy);
  %dy = [NaN dy];

  [~, indx] = max(ddy);
  indx = indx + 1;

  center = x(indx);

  [peak, peak_params] = get_peak(x - center, ddy);

  center = center + peak_params(1);
  peak_params(2) = 1.5*peak_params(2);

  dx = range(x);
  
  if (peak_params > dx/4)
    x_neg = (x < center - dx/4);
    x_pos = (x > center + dx/4);
  else
    x_neg = (x < center - peak_params(2));
    x_pos = (x > center + peak_params(2));
  end

  if (sum(x_neg) < 3)
    x_neg  = (x < center);
  end
  if (sum(x_pos) < 3)
    x_pos  = (x > center);
  end

  %val_all = [[zeros(length(x_neg), 2) ones(size(x_neg))]; ...
  %        [(center - x_pos(2:end)) ones(size(x_pos(2:end))) zeros(size(x_pos(2:end)))]] \ log(dy(2:end)).';

  %keyboard

  if (is_invag)
  
    vals3 = [[zeros(size(x(x_neg).')) (center - x(x_neg).') ones(size(x(x_neg))).']; ...
             [(center - x(x_pos).') zeros(size(x(x_pos).')) ones(size(x(x_pos).'))]] \ [(y(x_neg)).'; y(x_pos).'];

    params = real([center peak_params(2)/2 vals3.']);
  else

    valids = (dy < 0) & x_pos;

    if (sum(valids) < 3)
      valids = x_pos;
    end

    vals = [(center - x(valids).') ones(size(x(valids))).'] \ log(-dy(valids)).';
    %vals(2) = exp(vals(2)) / vals(1);
    %vals = real(vals);
    %vals2 = [(center - x(x_neg).') ones(size(x(x_neg))).'] \ (y(x_neg)).';

    vals3 = [[zeros(size(x(x_neg).')) (center - x(x_neg).') ones(size(x(x_neg))).']; ...
             [-(1 - exp(vals(1)*(center - x(x_pos).'))) zeros(size(x(x_pos).')) ones(size(x(x_pos).'))]] \ [(y(x_neg)).'; y(x_pos).'];

    params = real([center peak_params(2)/2 vals(1) vals3.']);
  end

  %figure;plot(x, dy);
  %figure;plot(x, y);
  %hold on;
  %plot(x, piecewise_background(params,x), 'r');
  
  %keyboard

  return;
end

function y = piecewise_background(p, x)
  
  if (size(p, 2) == 5)
    p = [p(:, [1:2]) NaN(size(p,1), 1) p(:, 3:end)];
  end

  problems = (any(p(:, 2:3) <= 0, 2) | any(~isfinite(p(:,[1:2 4:end])), 2));
  invags = (isnan(p(:, 3)) & ~problems);

  y = zeros(size(p, 1), length(x));

  p(:, 2) = 1.75*p(:,2);

  neg = bsxfun(@lt, x, p(:, 1) - p(:, 2));
  pos = bsxfun(@gt, x, p(:, 1) + p(:, 2));
  middle = (~neg & ~pos);

  knots = [-p(:, 5).*(-p(:, 2)) + p(:, 6), ...
           p(:, 6) - p(:, 4).*(1 - exp(-p(:, 3).*p(:, 2)))];

  tmp = bsxfun(@minus, p(:, 6), bsxfun(@times, p(:, 4), (1 - exp(bsxfun(@times, -p(:, 3), bsxfun(@minus, x, p(:, 1)))))));

  if (any(invags))
    knots(invags, 2) = -p(invags, 4).*p(invags,2) + p(invags,6);
    tmp(invags, :) = bsxfun(@plus, bsxfun(@times, -p(invags, 4), bsxfun(@minus, x, p(invags, 1))), p(invags, 6));
  end

  y(pos) = tmp(pos);

  tmp = bsxfun(@plus, bsxfun(@times, bsxfun(@rdivide, diff(knots, [], 2), 2*p(:, 2)), bsxfun(@minus, x, p(:, 1))), mean(knots, 2));
  y(middle) = tmp(middle);

  tmp = bsxfun(@plus, bsxfun(@times, -p(:, 5), bsxfun(@minus, x, p(:, 1))), p(:, 6));
  y(neg) = tmp(neg);

  y(problems, :) = 0;

  %y = p(4) * exp(-p(3) * (pos - p(1)).^2) + ...
  %    diff(knots)*middle + mean(knots) + ...
  %    -p(5) * neg + p(6);

  %y = p(3)*exp(-p(2)*(positive_value(x - p(1)))) + p(4) * positive_value(p(1) - x) + p(5);

  return;
end
