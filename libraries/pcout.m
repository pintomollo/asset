function outliers = pcout(pts, s)
% P. Filzmoser, R. Maronna, and M. Werner. Outlier identification in high dimensions.
% Computational Statistics and Data Analysis, Vol. 52, pp. 1694-1711, 2008

  if (nargin == 1)
    s = 0.25;
  end
  npts = size(pts, 1);

  rob_mean = nanmedian(pts);
  rob_var = mad(pts, 1);

  x = bsxfun(@rdivide, bsxfun(@minus, pts, rob_mean), rob_var);
  cov_mat = nancov(x);
  [eig_vect, eig_vals] = eig(cov_mat);

  eig_vals = diag(eig_vals);
  [eig_vals, indx] = sort(eig_vals, 'descend');
  eig_vect = eig_vect(:, indx);

  tot_var = sum(eig_vals);
  frac_exp = cumsum(eig_vals) / tot_var;
  last = find(frac_exp > 0.99, 1);
  eig_vect = eig_vect(:, 1:last);

  z = x*eig_vect;
  rob_mean = nanmedian(z);
  rob_var = mad(z, 1);
  z = bsxfun(@rdivide, bsxfun(@minus, z, rob_mean), rob_var);

  rob_mean = nanmedian(z);
  rob_var = mad(z, 1);
  rob_cov = nancov(z);

  centered = bsxfun(@minus, z, rob_mean);

  kurt = abs(nanmean(bsxfun(@rdivide, centered.^4, rob_var.^4)) - 3);
  kurt = kurt / sum(kurt);

  kcov = rob_cov .* (kurt' * kurt);
  kcov = pinv(kcov);

  quantiles = sqrt(chi2inv([0.25 0.5 0.99], last));

  mahala = NaN(npts, 1);
  for i=1:npts
    mahala(i) = sqrt(centered(i,:) * kcov * centered(i,:)');
  end
  mahala = mahala * (quantiles(2) / nanmedian(mahala));

  M = prctile(mahala, 100/3);
  c = nanmedian(mahala) + 2.5*mad(mahala, 1);
  w1 = biweight(mahala, c, M);

  rcov = pinv(rob_cov);
  for i=1:npts
    mahala(i) = sqrt(centered(i,:) * rcov * centered(i,:)');
  end
  mahala = mahala * (quantiles(2) / nanmedian(mahala));

  w2 = biweight(mahala, quantiles(3), quantiles(1));

  w = (w1 + s) .* (w2 + s) / ((1 + s)^2);
  outliers = (w < 0.25);

  return;
end

function w = biweight(d, c, M)

  w = zeros(size(d));
  w(d <= M) = 1;

  middle = (M < d) & (d < c);
  w(middle) = (1- ((d(middle) - M)./(c-M)).^2).^2;

  return;
end
