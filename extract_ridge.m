function values = extract_ridge(params, pos, dperp, rescale, opts)

  if (isstruct(params))
    mymovie = params;
    opts = pos;

    [nframes] = length(mymovie.data.quantification);

    for i = 1:nframes
      nimg = i;
      %nimg = 14 + i

      cortex = mymovie.data.quantification(nimg).carth;
      egg_dist = min(sqrt(bsxfun(@minus, cortex(:,1), mymovie.data.eggshell(nimg).carth(:,1).').^2 + bsxfun(@minus, cortex(:,2), mymovie.data.eggshell(nimg).carth(:,2).').^2), [], 2);
      egg_props = mymovie.data.quantification(nimg).eggshell;
      if (~any(isnan(egg_props)))
        autofluo = egg_props(1)*exp(-egg_dist.^2 ./ (2*egg_props(2)^2));
      else
        autofluo = zeros(size(egg_dist));
      end
      %if (opts.quantification.use_ruffles)
      %  [ cortex, rescale] = insert_ruffles(cortex, mymovie.data.ruffles(nimg).paths);
      %else
%        rescale = false(size(cortex, 1), 1);
      %end

      [dperp, dpos] = perpendicular_sampling(cortex, opts);

      if (size(mymovie.data.quantification(nimg).front, 1) ~= size(cortex, 1))
        error('Quantification and path do not correspond, you need to recompute both.')
      end

      rescale = mymovie.data.quantification(nimg).ruffles;

      mymovie.data.quantification(nimg).cortex = extract_ridge(mymovie.data.quantification(nimg).front, cortex, dperp, rescale, opts);
      mymovie.data.quantification(nimg).autofluorescence = autofluo;

      %pts = mymovie.data.quantification(nimg).front;

      %img = load_data(mymovie.data, nimg);
      %figure;imagesc(img);
      %size_img = size(img);
      %data = zeros(size_img);

      %centers = bsxfun(@times, pts(:,1), dperp) + cortex;

      %for j=1:size(pts, 1)
      %  j
      %  data = max(data, GaussMask2D(pts(j,2), size_img, centers(j,[2 1]), 0, 1)*pts(j,3));
      %end

      %figure;imagesc(data);
      %keyboard
    end

    values = mymovie;

    return;
  end

  if (numel(rescale) ~= size(params, 1))
    warning('Size of the ruffles does not match the size of the cortex, recompute !');

    return;
  end

  goods = (~any(isnan(params), 2));
  full_values = NaN(numel(goods), 1);

  params = params(goods, :);
  dperp = dperp(goods, :);
  pos = pos(goods, :);
  rescale = rescale(goods, :);

  npts = size(params, 1);
  values = NaN(npts,1);

  window_width = (opts.quantification.window_params / opts.pixel_size);
  gaussian_var = 1 / (2*(window_width^2));
  window_size = (1.79 * window_width)^2;

  %{
  dist = bsxfun(@minus, pos(:,1), pos(:,1).').^2 + bsxfun(@minus, pos(:,2), pos(:,2).').^2;
  dist(dist>window_size) = Inf;
  dist = exp(-dist * gaussian_var);
  dist = bsxfun(@rdivide, dist, sum(dist, 2));

  signal1 = sum(bsxfun(@times, dist, params(:,1).*exp(-params(:,2).^2 ./ (2*params(:,3).^2))),1).';

  centers = bsxfun(@times, params(:,2), dperp) + pos;

  dist = bsxfun(@minus, centers(:,1), pos(:,1).').^2 + bsxfun(@minus, centers(:,2), pos(:,2).').^2;
  dist(dist>window_size) = Inf;
  dist = exp(-dist * gaussian_var);
  dist = bsxfun(@rdivide, dist, sum(dist, 2));

  signal2 = sum(bsxfun(@times, dist, params(:,1)),1).';
  %}

  [lin_pos, tot] = carth2linear(pos);
  lin_pos = lin_pos*tot;

  dist = sqrt(bsxfun(@minus, lin_pos, lin_pos.').^2);
  dist = min(dist, abs(dist-tot)).^2;
  dist(dist>window_size) = Inf;
  dist = exp(-dist * gaussian_var);
  dist = bsxfun(@rdivide, dist, sum(dist, 2));

  %values = sum(bsxfun(@times, dist, params(:,1).*exp(-params(:,2).^2 ./ (2*params(:,3).^2))),1).';
  values = sum(bsxfun(@times, dist, params(:,1)),1).';

  %{

  centers = bsxfun(@times, params(:,2), dperp) + pos;
  window_width = (opts.quantification.window_params / opts.pixel_size);

  gaussian_var = 2*(window_width^2);
  mdist = median(sqrt(sum(diff(centers).^2, 2)));
  window_size = ceil(1.79 * window_width / mdist);

  mvar = median(params(:,3));

  centers = [centers(end-window_size+1:end, :); centers; centers(1:window_size,:)];
  params = [params(end-window_size+1:end, :); params; params(1:window_size,:)];

  indexes = [0:2*window_size];

  %factor = sqrt(2*pi);

  for i=1:npts
    window = centers(indexes+i, :);
    window_params = params(indexes+i, :);
    dist = exp(-(sum(bsxfun(@minus, window, pos(i,:)).^2, 2)) ./ gaussian_var);
    signal = exp(-(sum(bsxfun(@minus, window, pos(i,:)).^2, 2)) ./ (2*(window_params(:,3).^2))) .* window_params(:,1);
    
    % Maximum value
    %values(i) = max(signal);

    % Gaussian average
    values(i) = sum(signal .* dist / sum(dist));

    %%% MOLLINATION !!
    % Integral
    %values(i) = sum((window_params(:, 3) .* window_params(:,2)) .* (dist / sum(dist))) / mvar;
    % Amplitude
    %values(i) = sum(window_params(:, 3) .* (dist / sum(dist)));
  end
  %subplot(221)
  %scatter(centers(:,1), centers(:,2), 'b');
  %hold on
  %scatter(pos(:,1), pos(:,2), 'r');
  figure;
  plot([params(window_size+1:end-window_size, 1);params(window_size+1:end-window_size, 1)], 'b');
  hold on
  plot([values(1:npts, 1);values(1:npts, 1)], 'r');

  %}

  %figure;
  %plot(params(:,1), 'b');
  %hold on
  %plot(params(:,1).*exp(-params(:,2).^2 ./ (2*params(:,3).^2)), 'k');
  %plot(signal1, 'r');
  %plot(signal2, 'm');
  %plot(signal3, 'g');

  %keyboard

  %pause
  %values = signal1;

  [pos, len, val] = boolean_domains(rescale);
  ndom = length(pos);

  for i=1:ndom
    if (val(i) && len(i)>1)
      prev_indx = (i-1);
      if (prev_indx < 1)
        prev_indx = prev_indx + ndom;
      end
      next_indx = (i+1);
      if (next_indx > ndom)
        next_indx = next_indx - ndom;
      end

      if (len(prev_indx) > len(i))
        tmp_len = len(i);
      else
        tmp_len = len(prev_indx);
      end

      % need circular pixel picking
      left_vals = values(pos(prev_indx)+len(prev_indx) - [1:tmp_len]);

      if (len(next_indx) > len(i))
        tmp_len = len(i);
      else
        tmp_len = len(next_indx);
      end

      % need circular pixel picking
      right_vals = values(pos(next_indx) - 1 + [1:tmp_len]);

      if (ttest2(left_vals, right_vals))
        curr_vals = values([pos(i):pos(i)+len(i)-1]);
        half = round(len(i)/2);

        [lh, lp] = (ttest2(curr_vals, left_vals, 0.05, 'right'));
        [rh, rp] = (ttest2(curr_vals, right_vals, 0.05, 'right'));

        if (lh & rh)
          %if (lp > rp)
          if (ttest2(left_vals, right_vals, 0.05, 'right'));
            values(pos(i)+len(i) - [1:half]) = NaN;
          else
            values([pos(i):pos(i)+half-1]) = NaN;
          end
        elseif (lh)
          %values([pos(i):pos(i)+half-1]) = poissrnd(median(left_vals), 1, half);
          %values([pos(i):pos(i)+half-1]) = median(left_vals) + std(left_vals)*randn(1, half)/2;
          values([pos(i):pos(i)+half-1]) = NaN;
        elseif (rh)
          %values(pos(i)+len(i) - [1:half]) = poissrnd(median(right_vals), 1, half);
          %values(pos(i)+len(i) - [1:half]) = median(right_vals) + std(right_vals)*randn(1, half)/2;
          values(pos(i)+len(i) - [1:half]) = NaN;
        end

        %values(values < 0) = 0;
      end
    end
  end

  full_values(goods) = values;
  values = full_values;

  return;
end
