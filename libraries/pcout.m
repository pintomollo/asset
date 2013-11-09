function [outliers, w1, w2] = pcout(pts, w1coef, w2coef)
% Based on the R code
% P. Filzmoser, R. Maronna, and M. Werner. Outlier identification in high dimensions.
% Computational Statistics and Data Analysis, Vol. 52, pp. 1694-1711, 2008

  cs = 0.25;
  if (nargin < 3)
    w1coef = [0.33 2.5];
    w2coef = [0.25 0.99];
  else
    if (numel(w1coef) < 2)
      w1coef = [0.33 2.5];
    end
    if (numel(w2coef) < 2)
      w2coef = [0.25 0.99];
    end
  end
  npts = size(pts, 1);

  % PHASE 1:
  % Step 1: robustly sphere the data:
  rob_mean = nanmedian(pts);
  rob_var = 1.4826*mad(pts, 1);
  x = bsxfun(@rdivide, bsxfun(@minus, pts, rob_mean), rob_var);

  % Step 2: PC decomposition; compute p*, robustly sphere:
  [u, s, v] = svd(bsxfun(@minus, x, nanmean(x)));
  var_expl = (diag(s).^2) / (npts-1);
  frac_exp = cumsum(var_expl) / sum(var_expl);
  last = find(frac_exp > 0.99, 1);
  quantiles = sqrt(chi2inv([0.5 w2coef], last));

  z = x*v(:,1:last);
  rob_mean = nanmedian(z);
  rob_var = 1.4826*mad(z, 1);
  z = bsxfun(@rdivide, bsxfun(@minus, z, rob_mean), rob_var);

  % Step 3: compute robust kurtosis weights, transform to distances:
  kurt = abs(nanmean(z.^4) - 3);
  weighted = z * diag(kurt / sum(kurt));

  mahala = sqrt(sum(weighted.^2, 2));
  mahala = mahala * quantiles(1) / median(mahala);

  % Step 4: determine weights according to translated biweight:
  M = prctile(mahala, w1coef(1)*100);
  c = nanmedian(mahala) + w1coef(2)*1.4826*mad(mahala, 1);
  w1 = biweight(mahala, c, M);

  % PHASE 2:
  % Step 5: compute Euclidean norms of PCs and their distances:
  mahala = sqrt(sum(z.^2, 2));
  mahala = mahala * quantiles(1) / median(mahala);

  % Step 6: determine weight according to translated biweight:
  w2 = biweight(mahala, quantiles(3), quantiles(2));

  % Combine PHASE1 and PHASE 2: compute final weights:
  w = (w1 + cs) .* (w2 + cs) / ((1 + cs)^2);
  outliers = (w < 0.25);

  return;
end

function w = biweight(d, c, M)

  w = (1- ((d - M)./(c-M)).^2).^2;

  w(d > c) = 0;
  w(d <= M) = 1;

  return;
end
