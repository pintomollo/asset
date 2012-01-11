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

  if (opts.quantification.use_ruffles)
    if (opts.recompute | ~isfield(mymovie.(type), 'ruffles') | isempty(mymovie.(type).ruffles))
    %if (~isfield(mymovie.(type), 'ruffles') | isempty(mymovie.(type).ruffles))
      mymovie = find_ruffles(mymovie, opts);
      mymovie = follow_invaginations(mymovie, opts);
    elseif (~isfield(mymovie.(type).ruffles, 'paths') | isempty(mymovie.(type).ruffles(1).paths))
      do_it = true;
      for j=1:nframes
        if (~isempty(mymovie.(type).ruffles(j).paths))
          do_it = false;

          break;
        end
      end
      if (do_it)
        mymovie = follow_invaginations(mymovie, opts);
      end
    end
  end

  if (opts.recompute | ~isfield(mymovie.data, 'eggshell') | isempty(mymovie.data.eggshell))
    mymovie = duplicate_segmentation(mymovie, 'data', opts);
  end

  % Progress bar
  if (opts.verbosity > 0)
    hwait = waitbar(0,'Quantifying signal','Name','ASSET');
  end

  optims = optimset('Display', 'off');
      bin_dist = 64;
      bin_step = 1;
  
  for i=1:nframes

  bounds = [-Inf bin_step 0 0 -Inf 0 0; ...
             Inf Inf Inf Inf Inf Inf Inf];

  bounds_invag = [-Inf bin_step -Inf -Inf 0 0; ...
                   Inf Inf Inf Inf Inf Inf];


    nimg = i;
    %nimg = randi(nframes)
    %nimg = i + 5
    %nimg = 6

    cortex = mymovie.data.cortex(nimg).carth;
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

  %keyboard

      %ell_cortex = carth2elliptic(cortex, mymovie.data.centers(:, nimg), mymovie.data.axes_length(:, nimg), mymovie.data.orientations(1, nimg));

      %ell_pos = ell_cortex(:,1);
      %ell_pos = unique([ell_pos; 0; pi]);
      %ell_pos = ell_pos([diff([ell_pos; ell_pos(1)+2*pi]) > 1e-10]);

      %ell_cortex = interp_elliptic(ell_cortex, ell_pos);

      %cortex = elliptic2carth(ell_cortex, mymovie.data.centers(:, nimg), mymovie.data.axes_length(:, nimg), mymovie.data.orientations(1, nimg));

      %if (opts.quantification.use_ruffles)
      %  [ cortex, rescale] = insert_ruffles(cortex, mymovie.(type).ruffles(nimg).paths);
      %else
      %  rescale = false(size(cortex, 1), 1);
      %end

      %tmp_cortex = [cortex(end,:); cortex; cortex(1,:)];
      %dpos = (tmp_cortex(3:end, :) - tmp_cortex(1:end-2, :)) / 2;
      %dpts = diff(cortex(1:end-1, :));
%
      %dperp = [-dpos(:,2) dpos(:,1)];
      %nulls = all(dperp == 0, 2);
      %dperp(nulls, :) = dpts(nulls, :);

      %dperp = bsxfun(@rdivide, dperp, hypot(dperp(:,1), dperp(:,2))); 

      %dpos = [-bin_dist:bin_step:bin_dist];
      %nbins = length(dpos);

      %all_pos_x = bsxfun(@plus, dperp(:,1) * dpos, cortex(:,1));
      %all_pos_y = bsxfun(@plus, dperp(:,2) * dpos, cortex(:,2));
%
      %all_pos_x = all_pos_x(:);
      %all_pos_y = all_pos_y(:);

      %values = bilinear(img, all_pos_x, all_pos_y);
      %values = reshape(values, [size(cortex,1) nbins]);
      
      %%%%%%%%% BE CAREFULL OF MEMORY IF TOO LARGE !! 

      [values, dperp, dpos] = perpendicular_sampling(img, cortex, opts);
      [ph_values] = perpendicular_sampling(ph, cortex, dperp, dpos, opts);

      npts = size(values, 1);
      slopes = NaN(npts, 2);

      [params, bounds] = estimate_mean(dpos, ph_values(~rescale, :), false, bounds, optims);

      [params, junk, slopes(~rescale, 1)] = estimate_mean(dpos, values(~rescale, :), false, bounds, optims);
      bounds_invag(:, 1) = bounds(:, 1) + bounds(:, 2);
      bounds_invag(:, 2) = 2*bounds(:,2);
      [params_invag, bounds_invag, slopes(rescale, :)] = estimate_mean(dpos, values(rescale, :), true, bounds_invag, optims);

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
        %line_params = lsqcurvefit(@front, line_params, valid_dpos.', valid_line.', curr_bounds(1,:), curr_bounds(2, :), optims);
        line_params = fit_front(@front, line_params, valid_dpos.', valid_line.', curr_bounds(1,:), curr_bounds(2, :), optims);

        if (rescale(j))
          line_params = [line_params(1:2) NaN line_params(3:end)];
        end

        signal(j, :) = line_params;
      end

      mymovie.data.quantification(nimg).front = signal(:, [1 2 end]);
      mymovie.data.quantification(nimg).bkg = signal(:, [3:end-1]);
      mymovie.data.quantification(nimg).carth = cortex;

    else

      %tmp_cortex = [cortex(end,:); cortex; cortex(1,:)];
      %dpos = (tmp_cortex(3:end, :) - tmp_cortex(1:end-2, :)) / 2;

      %dperp = [-dpos(:,2) dpos(:,1)];
      %dperp = bsxfun(@rdivide, dperp, hypot(dperp(:,1), dperp(:,2))); 

      %signal = mymovie.data.quantification(nimg).front;

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
  %y = sum(imf(2:end, :), 1);
  if (size(imf, 1) > 2)
    smoothed = sum(imf(end-1:end, :), 1);
  else
    smoothed = imf(end, :);
  end

  %params = estimate_front(x, y, [], is_invag);
  tmp_params = [];
  if (isfinite(bounds(2,2)))
    tmp_params = mean(bounds(:,1:2), 1);
  end

  params = estimate_front(x, smoothed, tmp_params, is_invag);
  params = max(params, bounds(1,:));
  params = min(params, bounds(2,:));
  %params = lsqcurvefit(@front, params, x, y, bounds(1,:), bounds(2, :), optims);
  params = fit_front(@front, params, x, y, bounds(1,:), bounds(2, :), optims);

  if (nargout == 3)
    y = -mymean(diff(ys(:, x < params(1) - 3*params(2)), [], 2), 2);

    if (isempty(y))
      slopes = ones(size(ys, 1), 1) * params(end-1);
    else
      imf = emdc([], y);

      if (size(imf, 1) > 2)
        slopes = sum(imf(end-1:end, :), 1).';
      else
        slopes = imf(end, :).';
      end
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
  else
    slopes = [];
  end

  bounds(1, 1) = params(1) - 1.5*params(2);
  bounds(2, 1) = params(1) + 1.5*params(2);
  bounds(2, 1) = max(params(2) / 3, bounds(2,1));
  bounds(2, 2) = 2*params(2);

  return;
end

function new_params = estimate_front(x, y, params, is_invag)

  params_sigmoid = estimate_piecewise(x, y, params, is_invag);

  [estim, junk] = front([params_sigmoid 0], x);
  peak_range = get_peak(x, y-estim);
  params_gaussian = estimate_gaussian(x(peak_range), y(peak_range)-estim(peak_range));

  %new_params = [params_sigmoid(1) params_gaussian(2) params_sigmoid(3:end-1) params_sigmoid(end)+params_gaussian(end) params_gaussian(3:end-1)];
  new_params = [mean([params_sigmoid(1) params_gaussian(1)]) params_gaussian(2) params_sigmoid(3:end-1) params_sigmoid(end) params_gaussian(3:end-2) diff(params_gaussian([end end-1]))];

  return;
end

function [range, center] = get_peak(x, y)

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

  x = x - center;

  mins = find(dy(1:end-1) < 0 & dy(2:end) >= 0)+1;
  %mins = find(y(1:end-2) > y(2:end-1) & y(2:end-1) <= y(3:end)) + 1;
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

  if (isempty(best_min))
    best_min = 1;
  end

  range = mean([-x(lmin(best_min)), x(rmin(all_indx(best_min)))]);

  % Approx 1/20
  center = [center range/2.45];

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
    gauss = (0.9*exp(-((x-params(1)).^2) / (8*(params(2)^2))) + 0.1);
    ddy = -differentiator(x, dy .* gauss, true);
    dddy = differentiator(x, ddy, true);
    
    maxs = find(dddy(1:end-1) > 0 & dddy(2:end) <= 0) + 1;

    [~, indx] = min(abs(x(maxs) - params(1)));
    indx = maxs(indx);

    %if (isempty(indx))
   %   [~, indx] = max(ddy);
    %end
  else
    ddy = -differentiator(x, dy, true);
    [~, indx] = max(ddy);
  end

  if (isempty(indx))
    center = 0;
  else
    center = x(indx);
  end

  [peak, peak_params] = get_peak(x - center, ddy);

  center = center + peak_params(1);

  % Correction for estimating 2nd derivative
  peak_params(2) = sqrt(3)*peak_params(2);
  
  % Approx 1/5
  peak_dist = 1.79*peak_params(2);

  dx = range(x);
  
  if (peak_dist > dx/4)
    x_neg = (x < center - dx/4);
    x_pos = (x > center + dx/4);
  else
    x_neg = (x < center - peak_dist);
    x_pos = (x > center + peak_dist);
  end

  if (sum(x_neg) < 3)
    x_neg  = (x < center);
  end
  if (sum(x_pos) < 3)
    x_pos  = (x > center);
  end

  if (is_invag)
  
    vals3 = [[zeros(size(x(x_neg).')) (center - x(x_neg).') ones(size(x(x_neg))).']; ...
             [(center - x(x_pos).') zeros(size(x(x_pos).')) ones(size(x(x_pos).'))]] \ [(y(x_neg)).'; y(x_pos).'];

    params = [center peak_params(2) vals3.'];
  else

    %valids = (dy < 0) & x_pos;

    %if (sum(valids) < 3)
    %  valids = x_pos;
    %end

    %vals = real([(center - x(valids).') ones(size(x(valids))).'] \ log(-dy(valids)).');

    vals = real([(center - x(x_pos).') ones(size(x(x_pos))).'] \ log(-dy(x_pos)).');

    if (vals(1) < 0)
      valids = (x > center);
      vals = real([(center - x(valids).') ones(size(x(valids))).'] \ log(-dy(valids)).');

      if (vals(1) < 0)
        vals(1) = 1e-5;
      end
    end

    vals3 = [[zeros(size(x(x_neg).')) (center - x(x_neg).') ones(size(x(x_neg))).']; ...
             [-(1 - exp(vals(1)*(center - x(x_pos).'))) zeros(size(x(x_pos).')) ones(size(x(x_pos).'))]] \ [(y(x_neg)).'; y(x_pos).'];

    params = [center peak_params(2) vals(1) vals3.'];
  end

  return;
end
