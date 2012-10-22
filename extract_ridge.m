function values = extract_ridge(params, pos, dperp, rescale, opts)

  if (isstruct(params))
    mymovie = params;
    opts = pos;

    [nframes] = length(mymovie.data.quantification);

    for i = 1:nframes
      nimg = i;
      %nimg = 120 + i

      cortex = mymovie.data.quantification(nimg).carth;

      [dperp, dpos] = perpendicular_sampling(cortex, opts);

      if (size(mymovie.data.quantification(nimg).front, 1) ~= size(cortex, 1))
        error('Quantification and path do not correspond, you need to recompute both.')
      end

      rescale = isnan(mymovie.data.quantification(nimg).bkg(:,1));

      mymovie.data.quantification(nimg).cortex = extract_ridge(mymovie.data.quantification(nimg).front, cortex, dperp, rescale, opts);
    end

    values = mymovie;

    return;
  end

  goods = (~all(isnan(params), 2));
  params = params(goods, :);
  dperp = dperp(goods, :);
  pos = pos(goods, :);

  npts = size(params, 1);
  values = NaN(npts,1);
  safety = 1.5;

  centers = bsxfun(@times, params(:,1), dperp) + pos;

  gaussian_var = 2*(opts.quantification.window_params / opts.pixel_size);
  mdist = median(sqrt(sum(diff(centers).^2, 2)));
  window_size = ceil(safety * gaussian_var / mdist);

  mvar = median(params(:,2));

  centers = [centers(end-window_size+1:end, :); centers; centers(1:window_size,:)];
  params = [params(end-window_size+1:end, :); params; params(1:window_size,:)];

  indexes = [0:2*window_size];

  %factor = sqrt(2*pi);

  %% Might want to normalize the integral back to an amplitude ?!?

  for i=1:npts
    window = centers(indexes+i, :);
    window_params = params(indexes+i, :);
    dist = exp(-(sum(bsxfun(@minus, window, pos(i,:)).^2, 2)) ./ gaussian_var);

    %%% MOLLINATION !!

    % Integral
    %values(i) = sum((window_params(:, 3) .* window_params(:,2) .* factor) .* (dist / sum(dist))) / ;
    values(i) = sum((window_params(:, 3) .* window_params(:,2)) .* (dist / sum(dist))) / mvar;
    % Amplitude
    %values(i) = sum(window_params(:, 3) .* (dist / sum(dist)));
  end
  %subplot(221)
  %scatter(centers(:,1), centers(:,2), 'b');
  %hold on
  %scatter(pos(:,1), pos(:,2), 'r');
  %subplot(222)
  %plot(params(:,3));
  %subplot(224)
  %plot(values);

  %pause

  [pos, len, val] = boolean_domains(rescale);
  ndom = length(pos);

  for i=1:ndom
    if (val(i))
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

  return;
end
