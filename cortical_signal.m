function mymovie = cortical_signal(mymovie, opts)

  type = opts.segmentation_type;
  if (strncmp(type, 'markers', 7))
    [nframes imgsize] = size_data(mymovie.cortex);
  else
    [nframes imgsize] = size_data(mymovie.(type));
  end

  if (~opts.recompute && isfield(mymovie.data, 'quantification') && ...
     ~isemtpy(mymovie.data.quantification) && length(mymovie.data.quantification) == nframes)
    
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
  
  for i=1:nframes
    %nimg = randi(nframes)
    nimg = i;
    %nimg = i + 73
    %nimg = 52
    %keyboard

    cortex = mymovie.data.cortex(nimg).carth;

    if (opts.recompute|(~isfield(mymovie.data, 'quantification'))|(length(mymovie.data.quantification) < nimg)|~isfield(mymovie.data.quantification(nimg), 'front')|isempty(mymovie.data.quantification(nimg).front))
      img = double(load_data(mymovie.data, nimg));

      ell_cortex = carth2elliptic(cortex, mymovie.data.centers(:, nimg), mymovie.data.axes_length(:, nimg), mymovie.data.orientations(1, nimg));

      ell_pos = ell_cortex(:,1);
      ell_pos = unique([ell_pos; 0; pi]);
      ell_pos = ell_pos([diff([ell_pos; ell_pos(1)+2*pi]) > 1e-10]);

      ell_cortex = interp_elliptic(ell_cortex, ell_pos);

      cortex = elliptic2carth(ell_cortex, mymovie.data.centers(:, nimg), mymovie.data.axes_length(:, nimg), mymovie.data.orientations(1, nimg));

      if (opts.quantification.use_ruffles)
        [ cortex, rescale] = insert_ruffles(cortex, mymovie.(type).ruffles(nimg).carth, mymovie.(type).ruffles(nimg).paths);
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

      dpos = [-nbins:0.5:nbins];
      nbins = length(dpos);
      sampling_index = [1:nbins];

      all_pos_x = bsxfun(@plus, dperp(:,1) * dpos, cortex(:,1));
      all_pos_y = bsxfun(@plus, dperp(:,2) * dpos, cortex(:,2));

      all_pos_x = all_pos_x(:);
      all_pos_y = all_pos_y(:);

      values = bilinear(img, all_pos_x, all_pos_y);
      values = reshape(values, [size(cortex,1) nbins]);

      projection = mymean(values, 1);
      valids = ~isnan(projection);

      valid_projection = projection(valids);
      valid_dpos = dpos(valids);

      imf = emdc(sampling_index(valids), valid_projection);
      valid_projection = sum(imf(2:end, :), 1);

      params1 = estimate_sigmoid(valid_dpos, valid_projection);
      estim1 = sigmoid(params1, valid_dpos);

      peak_range = get_peak(valid_dpos, valid_projection-estim1);
      params2 = estimate_gaussian(valid_dpos(peak_range), valid_projection(peak_range)-estim1(peak_range));
      estim2 = gaussian([params2(1:3) 0], valid_dpos);

      params3 = estimate_sigmoid(valid_dpos, valid_projection - estim2);
      estim3 = sigmoid(params3, valid_dpos);

      peak_range = get_peak(valid_dpos, valid_projection - estim3);
      params4 = estimate_gaussian(valid_dpos(peak_range), valid_projection(peak_range)-estim3(peak_range));

      params5 = params3;
      params5(4) = params5(4) + params4(4);
      estim5 = sigmoid(params5, dpos);

      peak_range = correct_peak(peak_range, valids);

      npts = size(values, 1);

      signal = NaN(npts, 4); 
      smoothed = NaN(size(values));
      for j=1:npts
        line = values(j,:);
        valids = ~isnan(line);

        if (~any(valids))
          all_pos_x = reshape(all_pos_x, [size(cortex,1) nbins]);
          all_pos_y = reshape(all_pos_y, [size(cortex,1) nbins]);

          [all_pos_x(j, :); all_pos_y(j, :)]

          j

          continue;
        end

        valid_line = line(valids);

        valid_range = valids & peak_range;

        tmp_val = emdc(sampling_index(valids), valid_line);
        smoothed(j, valids) = sum(tmp_val(end-1:end, :), 1);
        signal(j, :) = estimate_gaussian(dpos(valid_range), smoothed(j, valid_range)-estim5(valid_range));

        if (rescale(j))
          signal(j, 3) = signal(j, 3) / 2;
        end
      end

      mymovie.data.quantification(nimg).front = signal;
      mymovie.data.quantification(nimg).bkg = [params5 params4(1:3)];
      mymovie.data.quantification(nimg).carth = cortex(:, :);

      if (false)
      estim4 = gaussian([params4(1:3) 0], dpos);
      clim = [min(values(:)) max(values(:))];
      figure;imagesc(gaussian(signal, dpos));
      figure;imagesc(bsxfun(@plus, gaussian(signal, dpos), estim5), clim);
      figure;imagesc(values, clim);
      figure;imagesc(smoothed);
      figure;plot(dpos, projection);hold on;plot(dpos(peak_range), projection(peak_range), 'r');
      plot(dpos, estim5, 'k');plot(dpos, estim4, 'g');
      plot(dpos, estim5 + estim4, 'm')
      %plot(dpos, estim3, 'm');plot(dpos, estim2, 'y');
      %plot(dpos, estim1, 'c');
      keyboard
      end
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

  if (false)
  %if (opts.recompute|~isfield(mymovie.data, 'domain')|isempty(mymovie.data.domain))
    img = gather_quantification(mymovie, opts);
    img = imnorm(img);
    opts.quantification.params.init = sub2ind(size(img,2)*[1 1], [1:size(img,2)],[1:size(img,2)]);
    path = dynamic_prog_2d(img, opts.quantification.params, @weight_domain_borders, opts.quantification.weights, opts);
    mymovie.data.domain = path / size(img, 2);
  end

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

function range = get_peak(x, y)

  % Get only the peak
  dy = diff(y);

  ymax = find(dy(1:end-1) > 0 & dy(2:end) <= 0)+1;
  if (isempty(ymax))
    center = 0;
  else
    [~, center] = min(abs(x(ymax)));
    center = x(ymax(center));
  end

  x = x - center;

  %lmin = find(x(2:end-1) < 0 & dy(1:end-1) < 0 & dy(2:end) >= 0)+1;
  %rmin = find(x(2:end-1) > 0 & dy(1:end-1) <= 0 & dy(2:end) > 0)+1;

  mins = find(dy(1:end-2) > dy(2:end-1) & dy(2:end-1) <= dy(3:end)) + 1;
  maxs = [];
  lmin = mins(x(mins) < 0);
  if (isempty(lmin))
    maxs = find(dy(1:end-2) < dy(2:end-1) & dy(2:end-1) >= dy(3:end)) + 1;
    lmin = maxs(x(maxs) < 0);
    if (isempty(lmin))
      lmin = 1;
    end
  end
  rmin = mins(x(mins) > 0);
  if (isempty(rmin))
    if (isempty(maxs))
      maxs = find(dy(1:end-2) < dy(2:end-1) & dy(2:end-1) >= dy(3:end)) + 1;
    end
    rmin = maxs(x(maxs) > 0);
    if (isempty(rmin))
      rmin = length(x);
    end
  end

  width = abs(bsxfun(@minus, x(rmin), x(lmin).'));
  dist = abs(bsxfun(@plus, x(rmin), x(lmin).'));
  dist(width < 30) = Inf;
  [dist, indexes] = min(dist, [], 2);
  best_min = find(dist(2:end) <= dist(1:end-1), 1, 'last') + 1;
  if (isempty(best_min))
    best_min = length(dist);
  end

  range = min(-x(lmin(best_min)), x(rmin(indexes(best_min))));
  range = (x >= -range & x <= range);

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

function y = sigmoid(p, x)

  y = bsxfun(@plus, bsxfun(@times, ((1 + erf(bsxfun(@rdivide, bsxfun(@plus, -x, p(:, 1)), p(:, 2)*sqrt(2)))) ./ 2), p(:, 3)), p(:, 4));

  return;
end

function y = gaussian(p, x)

  y = bsxfun(@plus, bsxfun(@times, p(:, 3), exp(bsxfun(@rdivide, -(bsxfun(@minus,x,p(:,1)).^2), (2*(p(:, 2).^2))))), p(:, 4));

  return;
end
