function mymovie = cortical_signal(mymovie, opts)

  nframes = size_data(mymovie.cortex);

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
  
  for i=1:nframes
    %nimg = randi(nframes)
    nimg = i;
    %nimg = 38;

    cortex = mymovie.data.cortex(nimg).carth;

    if (opts.recompute|(length(mymovie.data.quantification) < nimg)|~isfield(mymovie.data.quantification(nimg), 'front')|isempty(mymovie.data.quantification(nimg).front))
      img = double(load_data(mymovie.data, nimg));

      ell_cortex = abs(carth2elliptic(cortex, mymovie.data.centers(:, nimg), mymovie.data.axes_length(:, nimg), mymovie.data.orientations(1, nimg)));
      indx = find(ell_cortex(:,1)==min(ell_cortex(:,1)), 1);

      if (indx ~= 1)
        cortex = cortex([indx:end 1:indx-1], :);
        mymovie.data.cortex(nimg).carth = cortex;
      end

      tmp_cortex = [cortex(end,:); cortex; cortex(1,:)];
      dpos = (tmp_cortex(3:end, :) - tmp_cortex(1:end-2, :)) / 2;

      dperp = [-dpos(:,2) dpos(:,1)];
      dperp = bsxfun(@rdivide, dperp, hypot(dperp(:,1), dperp(:,2))); 

      nbins = 64;

      dpos = [-nbins:0.5:nbins];
      nbins = length(dpos);

      all_pos_x = bsxfun(@plus, dperp(:,1) * dpos, cortex(:,1));
      all_pos_y = bsxfun(@plus, dperp(:,2) * dpos, cortex(:,2));

      all_pos_x = all_pos_x(:);
      all_pos_y = all_pos_y(:);

      values = bilinear(img, all_pos_x, all_pos_y);
      values = reshape(values, [size(cortex,1) nbins]);

      projection = mean(values, 1);
      imf = emdc([], projection);
      if (any(isnan(imf(:))))
      %  imf = emd([], projection);
        'NaN proj !'
        keyboard
      end
      projection = imf(end, :);

      params1 = estimate_sigmoid(dpos, projection);
      estim1 = sigmoid(params1, dpos);

      peak_range = get_peak(dpos, projection-estim1);
      params2 = estimate_gaussian(dpos(peak_range), projection(peak_range)-estim1(peak_range));
      estim2 = gaussian([params2(1:3) 0], dpos);

      %params3 = nlinfit(dpos, projection - estim2, @sigmoid, params1);
      %estim3 = sigmoid(params3, dpos);

      %params4 = nlinfit(dpos, projection - estim3, @gaussian, params2);
      %estim4 = gaussian(params4, dpos);

      params3 = estimate_sigmoid(dpos, projection - estim2);
      estim3 = sigmoid(params3, dpos);

      peak_range = get_peak(dpos, projection - estim3);
      params4 = estimate_gaussian(dpos(peak_range), projection(peak_range)-estim3(peak_range));
      estim4 = gaussian([params4(1:3) 0], dpos);

      params5 = params3;
      params5(4) = params5(4) + params4(4);
      estim5 = sigmoid(params5, dpos);

      npts = size(values, 1);

      signal = NaN(npts, 4); 
      dpos_range = dpos(peak_range);
      estim_range = estim5(peak_range);
      %signal2 = NaN(npts, 4); 
      %signal3 = NaN(npts, 4); 
      %signal4 = NaN(npts, 4); 
      for j=1:npts
        %signal(j, :) = nlinfit(dpos, values(j, :) - estim3, @gaussian, params4);
        %signal2(j, :) = estimate_gaussian(dpos(peak_range), values(j, peak_range)-estim3(peak_range));
        line = values(j,:);
        if (any(isnan(line)))
          disp('?!?!?');
        end
        tmp_val = emdc([],line);
        %while (any(isnan(tmp_val(:))))
        %  tmp_val = emdc([], line);
        %end
        signal(j, :) = estimate_gaussian(dpos_range, tmp_val(end, peak_range)-estim_range);
        %signal4(j, :) = nlinfit(dpos, values(j, :) - estim3, @gaussian, signal2(j,:));
      end
      mymovie.data.quantification(nimg).front = signal;
      mymovie.data.quantification(nimg).bkg = [params5 params4(1:3)];

      %figure;imagesc(gaussian(signal, dpos));
      %keyboard
    else

      tmp_cortex = [cortex(end,:); cortex; cortex(1,:)];
      dpos = (tmp_cortex(3:end, :) - tmp_cortex(1:end-2, :)) / 2;

      dperp = [-dpos(:,2) dpos(:,1)];
      dperp = bsxfun(@rdivide, dperp, hypot(dperp(:,1), dperp(:,2))); 

      signal = mymovie.data.quantification(nimg).front;
    end

    %figure;plot(dpos,projection);
    
    if (opts.recompute|length(mymovie.data.quantification) < nimg|~isfield(mymovie.data.quantification(nimg), 'cortex')|isempty(mymovie.data.quantification(nimg).cortex))
      values = extract_ridge(signal, cortex, dperp, opts);
      mymovie.data.quantification(nimg).cortex = values;
    end

    %hold on;plot(dpos,estim1,'r');
    %plot(dpos, estim3, 'g')

    %figure;plot(dpos,projection-estim1);
    %hold on;plot(dpos,estim2,'r');
    %plot(dpos, estim4, 'g')

    %index = (dpos < -8 | dpos > 8);
    %clean_proj = projection;
    %gap_length = sum(~index);
    %clean_proj(~index) =  ([1:gap_length] / (gap_length + 1));
    %params2 = estimate_sigmoid(dpos(index), projection(index));

    %figure;imagesc(values);
    %figure;imagesc(gaussian(signal, dpos));

    %keyboard
  end

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

  %bkg = min(y);
  %ampl = max(y) - bkg;

  %y = (y - bkg) / ampl;
  %[~, center] = min(abs(y - 0.5));
  %center = x(center);
  %x = (x - center);

  %index = (y ~= 0 & y ~= 1 & x~= 0);

  %slope = mean(log((1 ./ y(index)) - 1) ./ x(index));

  %params = [center, slope, ampl, bkg];

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

  lmin = find(x(2:end-1) < 0 & dy(1:end-1) < 0 & dy(2:end) >= 0, 1, 'last')+1;
  if (isempty(lmin))
    lmin = 1;
  end
  rmin = find(x(2:end-1) > 0 & dy(1:end-1) <= 0 & dy(2:end) > 0, 1, 'first')+1;
  if (isempty(rmin))
    rmin = length(x);
  end

  range = min(-x(lmin), x(rmin));
  range = (x >= -range & x <= range);

  return;
end

function params = estimate_gaussian(x,y)

  if (any(isnan(y)))
    'NaNssss'
    beep;beep;
    keyboard
  end

  orig_y = y;
  y = y - min(y);
  
  y = y / sum(y);
  center = sum(x .* y);

  x = (x - center).^2;
  sigma = sum(y .* x);

  % Try to catch the possible warning (lastwarn)
  % and adapt to it, either reduce #params or try pinv
  % check 14-110311_2!!! 24-170211_1 !! 24-090311_4


  %ampl = mean(orig_y(range) ./ exp(-x / (2*sigma)));
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

  %y = bsxfun(@plus, bsxfun(@rdivide, p(:, 3), (1+exp(bsxfun(@times, p(:, 2), bsxfun(@minus, x, p(:, 1)))))), p(:, 4));
  y = bsxfun(@plus, bsxfun(@times, ((1 + erf(bsxfun(@rdivide, bsxfun(@plus, -x, p(:, 1)), p(:, 2)*sqrt(2)))) ./ 2), p(:, 3)), p(:, 4));

  return;
end

function y = gaussian(p, x)

  y = bsxfun(@plus, bsxfun(@times, p(:, 3), exp(bsxfun(@rdivide, -(bsxfun(@minus,x,p(:,1)).^2), (2*(p(:, 2).^2))))), p(:, 4));

  return;
end
