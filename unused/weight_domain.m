function weight = weight_domain(img, params)

  alpha = params.alpha;
  beta = params.beta;
  gamma = params.gamma;

  [nrows,npts] = size(img);
  nsize = floor(npts/2) + 1;
  
  img = gaussian_mex(img, 0.6);
  minmax = prctile(img(:), [0.1 99.9]);
  img = imnorm(img, minmax(1), minmax(2));

  weight = zeros(nrows, npts, nsize);

  derivatives = 0.5 * imnorm(differentiator(img, 2));

  left_indexes = [1:npts];
  right_indexes = left_indexes;

  valids = isfinite(img);
  img(~valids) = 0;
  dist = cumsum(valids, 2);
  max_dist = dist(:,end);

  dist(~valids) = NaN;
  vals = cumsum(img, 2);
  max_vals = vals(:, end);

  img(~valids) = NaN;

  sums = [zeros(nrows, 1) vals(:, 1:end-1)];
  dists = [zeros(nrows, 1) dist(:, 1:end-1)];

  for i=1:nsize

    if (i == 1)
      tmp_deriv = 4*abs(derivatives - 0.25);
      tmp_pixels = img;

      domain_value = repmat(1 - (max_vals ./ max_dist), 1, npts);
    else

      left_indexes = circshift(left_indexes, [0 1]);

      domain_size = dists(:, right_indexes) - dists(:, left_indexes);
      is_inverted = (domain_size < 0);
      if (any(is_inverted(:)))
        tmp_size = bsxfun(@plus, max_dist, domain_size);
        domain_size(is_inverted) = tmp_size(is_inverted);
      end

      domain_integral = sums(:, right_indexes) - sums(:, left_indexes);
      complement_integral = bsxfun(@plus, max_vals, domain_integral);
      domain_integral(is_inverted) = complement_integral(is_inverted);

      domain = domain_integral ./ domain_size;
      outside = bsxfun(@minus, max_vals, domain_integral) ./ bsxfun(@minus, max_dist, domain_size);

      domain_value = beta .* (1-domain) + (1-beta) .* outside;

      right_indexes = circshift(right_indexes, [0 -1]);

      tmp_deriv = 0.5 - derivatives(:, left_indexes) + derivatives(:, right_indexes);
      tmp_pixels = (img(:, left_indexes) + img(:, right_indexes)) ./ 2;
    end

    pixel_value = gamma .* tmp_deriv + (1-gamma) .* tmp_pixels;
    weight(:, :, i) = alpha .* domain_value + (1-alpha) .* pixel_value;
  end

  weight(isnan(weight)) = Inf;
  
  return;
end
