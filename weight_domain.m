function weight = weight_domain(img, params)

  alpha = params.alpha;
  beta = params.beta;
  gamma = 1/params.gamma;

  [nrows,npts] = size(img);
  nsize = floor(npts/2) + 1;
  img = imnorm(img, [], [], 'row');
  img = gaussian_mex(img, 0.6);

  all_derivatives = zeros(nrows, npts, nsize);
  domain_value = zeros(nrows, npts, nsize);

  left_deriv = zeros(nrows, npts);
  right_deriv = zeros(nrows, npts);
  left_sum = zeros(nrows, npts);
  right_sum = zeros(nrows, npts);
  left_size = zeros(nrows, npts);
  right_size = zeros(nrows, npts);

  %line = img(40, :);

  %valids = ~isnan(line);
  %line(~valids) = 0;
  %dist = cumsum(valids);

  %domain_size = bsxfun(@minus, dist, dist.') + 1;
  %vals = cumsum(line);
  %domain_integral = bsxfun(@minus, vals, [0 vals(1:end-1)].');
  
  %domain = domain_integral ./ domain_size;
  %outside = (vals(end) - domain_integral) ./ (npts - domain_size);
  %domain_vals = beta * (1 - domain) + (1-beta)*outside;

  %valids_deriv = (valids & valids([2:end 1]));
  %dist = cumsum(valids_deriv);

  %derivatives = diff([line line(1)]) ./ diff([0 dist]);
  %derivatives = 1 ./ (1 + exp(-delta * derivatives));
  
  %left_deriv = gamma * (1-derivatives);
  %right_deriv = (1-gamma) * derivatives;
  %derivatives = bsxfun(@plus, left_deriv.', right_deriv);

  %weight = alpha * domain_vals + (1-alpha) * derivatives;
  %weight(domain_size <= 0) = Inf;
  %weight(isnan(weight)) = Inf;

  %figure;imagesc(weight);

  mirror = img(:,[end-2:end 1:end 1:3]);
  derivatives = (39*(mirror(:,5:end-2) - mirror(:,3:npts+2)) + 12*(mirror(:,6:end-1) - mirror(:,2:npts+1)) -5*(mirror(:, 7:end) - mirror(:, 1:npts))) / 96;
  derivatives = 1 ./ (1 + exp(-gamma * derivatives));
  left_deriv = 0.5 * (1-derivatives);
  right_deriv = 0.5 * derivatives;

  %all_derivatives = bsxfun(@plus, left_deriv, reshape(right_deriv, [nrows, 1, npts]));

  valids = ~isnan(img);
  img(~valids) = 0;
  dist = cumsum(valids, 2);
  max_dist = dist(:,end);

  dist(~valids) = NaN;
  vals = cumsum(img, 2);
  max_vals = vals(:, end);

  left_sum = [zeros(nrows, 1) vals(:, 1:end-1)];
  right_sum = vals;
  left_dist = [zeros(nrows, 1) dist(:, 1:end-1)];
  right_dist = dist;

  for i=1:nsize


    if (i == 1)
      derivatives(derivatives < 0.5) = 1 - derivatives(derivatives < 0.5);
      all_derivatives(:, :, i) = derivatives;
    else
      all_derivatives(:, :, i) = left_deriv + right_deriv;
    end
    left_deriv = circshift(left_deriv, [0 1]);
    right_deriv = circshift(right_deriv, [0 -1]);

    domain_size = right_dist - left_dist;
    is_inverted = (domain_size < 0);

    if (any(is_inverted(:)))
      tmp_size = bsxfun(@plus, max_dist, domain_size);
      domain_size(is_inverted) = tmp_size(is_inverted);
    end

    domain_integral = right_sum - left_sum;
    complement_integral = bsxfun(@plus, max_vals, domain_integral);
    domain_integral(is_inverted) = complement_integral(is_inverted);

    domain = domain_integral ./ domain_size;
    outside = bsxfun(@minus, max_vals, domain_integral) ./ bsxfun(@minus, max_dist, domain_size);
    domain_value(:, :, i) = beta .* (1-domain) + (1-beta) .* outside;

    left_dist = circshift(left_dist, [0 1]);
    right_dist = circshift(right_dist, [0 -1]);

    left_sum = circshift(left_sum, [0 1]);
    right_sum = circshift(right_sum, [0 -1]);
  end

  %domain_size = bsxfun(@minus, reshape(dist, [nrows 1 npts]), dist) + 1;
  %is_inverted = (domain_size < 0);
  %domain_size(is_inverted) = domain_size(is_inverted) + npts;

  %vals = cumsum(img, 2);
  %domain_integral = bsxfun(@minus, reshape(vals, [nrows, 1, npts]), [zeros(nrows, 1) vals(:, 1:end-1)]);
  %complement_integral = bsxfun(@plus, vals(:, end), domain_integral);
  %domain_integral(is_inverted) = complement_integral(is_inverted);
  
  %domain = domain_integral ./ domain_size;
  %outside = bsxfun(@minus, vals(:, end), domain_integral) ./ (npts - domain_size);
  %domain_value = beta .* (1-domain) + (1-beta).*outside;
  
  %keyboard

  weight = (alpha * domain_value) + ((1-alpha) * all_derivatives);
  %weight(domain_size <= 0) = Inf;
  weight(isnan(weight)) = Inf;

  %figure;imagesc(squeeze(weight(40,:,:)));

  %keyboard
  
  return;
  
  
  domain_prop = domain_size ./ (npts + 1 - domain_size);

  vals = cumsum(line);
  domain_integral = bsxfun(@minus, vals, vals.');
  
  relative_mean = ((vals(end) - domain_integral) ./ domain_integral) .* domain_prop;
  relative_mean(domain_prop <= 0) = Inf;
  relative_mean(isnan(relative_mean)) = Inf;

  pts = [1:npts];
  pts = pts(valids);
  smoothed = emdc(pts, line(valids));
  smoothed = sum(smoothed(2:end, :));
  derivatives = diff([smoothed smoothed(1)]) ./ diff([pts pts(1)+npts]);


  keyboard

  %sums = cumsum(img, 2);

  return;
end
