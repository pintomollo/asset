function [all_sumvals, group] = mysum(all_vals, dim, indexes)

  all_sumvals = [];

  if (nargin == 2)
    indexes = ones(size(all_vals));
  elseif (nargin == 1)
    dim = 1;
    indexes = ones(size(all_vals));
  end

  if (isempty(dim))
    dim = 1;
  end
  groups = unique(indexes(:)).';

  for g = groups
    if (length(groups) == 1)
      vals = all_vals;
    else
      vals = all_vals(indexes == g);
    end

    nans = isnan(vals);
    nelems = sum(~nans, dim);
    vals(nans) = 0;

    sumvals = sum(vals, dim);
    sumvals(nelems == 0) = NaN;

    all_sumvals = cat(dim, all_sumvals, sumvals);

    if (nargout == 1)
      continue;
    end

    old_size = ones(1,ndims(vals));
    old_size(dim) = size(vals, dim);

  end

  return;
end
