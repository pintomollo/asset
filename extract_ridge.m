function values = extract_ridge(params, pos, dperp, rescale, opts)

  if (isstruct(params))
    mymovie = params;
    opts = pos;

    [nframes] = length(mymovie.data.quantification);

    for i = 1:nframes
      nimg = i;

      cortex = mymovie.data.cortex(nimg).carth;
      if (opts.quantification.use_ruffles)
        [ cortex, rescale] = insert_ruffles(cortex, mymovie.(type).ruffles(nimg).paths);
      else
        rescale = false(size(cortex, 1), 1);
      end

      [dperp, dpos] = perpendicular_sampling(cortex, opts);

      if (size(mymovie.data.quantification(nimg).front, 1) ~= size(cortex, 1))
        error('Quantification and path do not correspond, you need to recompute both.')
      end

      mymovie.data.quantification(nimg).cortex = extract_ridge(mymovie.data.quantification(nimg).front, cortex, dperp, rescale, opts);
    end

    values = mymovie;

    return;
  end

  window_size = 15;

  npts = size(params, 1);
  values = NaN(npts,1);

  centers = bsxfun(@times, params(:,1), dperp) + pos;

  centers = [centers(end-window_size+1:end, :); centers; centers(1:window_size,:)];
  params = [params(end-window_size+1:end, :); params; params(1:window_size,:)];

  indexes = [0:2*window_size];
  gaussian_var = opts.quantification.window_params / opts.pixel_size;
  gaussian_var = 2*(gaussian_var^2);

  factor = sqrt(2*pi);

  %% Might want to normalize the integral back to an amplitude ?!?

  for i=1:npts
    window = centers(indexes+i, :);
    window_params = params(indexes+i, :);
    dist = exp(-sqrt(sum(bsxfun(@minus, window, pos(i,:)).^2, 2)) ./ gaussian_var);

    % Integral
    values(i) = sum((window_params(:, 3) .* window_params(:,2) .* factor) .* (dist / sum(dist)));
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



  return;
end
