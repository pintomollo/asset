function [mu, sigma, is_pos] = folded_distribution(values, probs, dim)

  if (nargin == 1)
    probs = ones(size(values));
    dim = 1;
  elseif (nargin == 2)
    if (numel(probs) == 1)
      dim = probs;
      probs = ones(size(values));
    else
      dim = 1;
    end
  end

  sizes = size(values);
  perm_dim = [1:length(sizes)];
  perm_dim(dim) = 1;
  perm_dim(1) = dim;

  values = permute(values, perm_dim);
  values = reshape(values, sizes(dim), []);

  probs = permute(probs, perm_dim);
  probs = reshape(probs, sizes(dim), []);
  probs(~isfinite(probs)) = 0;

  total = sum(probs);
  probs(:, total == 0) = 1;
  probs = bsxfun(@rdivide, probs, sum(probs));

  m1 = bsxfun(@times, values, probs);
  m2 = bsxfun(@times, values.^2, probs);

  m1(~isfinite(m1)) = 0;
  m2(~isfinite(m2)) = 0;

  m1 = sum(m1);
  m2 = sum(m2);

  a = (pi - 2) / 2;
  b = m1*(2-pi);
  c = m1.^2 * (pi/2) - m2;
  
  mu1 = (-b + sqrt(b.^2 -4*a.*c)) ./ (2*a);
  sigma1 = (mu1 - m1)*sqrt(pi/2);

  mu2 = (-b - sqrt(b.^2 -4*a.*c)) ./ (2*a);
  sigma2 = (m1 - mu2)*sqrt(pi/2);

  p1 = sum(log(folded(values, mu1, sigma1, true) + 1e-10));
  p2 = sum(log(folded(values, mu2, sigma2, false) + 1e-10));
  is_pos = (p2 > p1);

  mu = mu1;
  mu(is_pos) = mu2(is_pos);
  sigma = sigma1;
  sigma(is_pos) = sigma2(is_pos);

  mu = ipermute(mu, perm_dim);
  sigma = ipermute(sigma, perm_dim);
  is_pos = ipermute(is_pos, perm_dim);

  return;
end

function y = folded(x, mu, sigma, greater)

  if (nargin == 3)
    greater = true;
  end

  %y = (2/(sigma * sqrt(2*pi))) .* exp(-(x-mu).^2/(2*sigma.^2));
  y = bsxfun(@times, (2./(sigma * sqrt(2*pi))), exp(bsxfun(@rdivide, -bsxfun(@minus, x, mu).^2, (2*sigma.^2))));

  if (greater)
    y(bsxfun(@gt, x, mu)) = 0;
  else
    y(bsxfun(@lt, x, mu)) = 0;
  end

  return;
end
