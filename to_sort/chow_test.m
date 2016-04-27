function p_all = chow_test(set1, set2)

  [ess, N] = get_residuals([set1;set2]);
  [ess1, N1] = get_residuals(set1);
  [ess2, N2] = get_residuals(set2);

  ndims = size(set1, 2) - 1;

  Chow = ((ess - (ess1 + ess2)) / ndims) / ((ess1 + ess2) / (N1 + N2 - 2*ndims));

  p_all = fcdf(Chow, ndims, N1 + N2 - 2*ndims);

  return;
end

function [ess, npts] = get_residuals(pts)

  x = pts(:,1:end-1);
  y = pts(:,end);

  [b, bint, r, rint] = regress(y, x);

  ess = sum(r.^2);
  npts = length(y);

  return;
end
