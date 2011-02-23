function [all_mvals, all_svals, groups] = mymean(all_vals, dim, indexes)

  all_mvals = [];
  all_svals = [];

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

    mvals = sum(vals, dim) ./ nelems;
    mvals(nelems == 0) = NaN;

    all_mvals = cat(dim, all_mvals, mvals);

    if (nargout == 1)
      continue;
    end

    old_size = ones(1,ndims(vals));
    old_size(dim) = size(vals, dim);

    svals = (vals - repmat(mvals, old_size)).^2;
    svals(nans) = 0;
    svals = sqrt(sum(svals, dim) ./ (nelems - 1));
    svals(nelems == 0) = NaN;
    svals(nelems == 1) = 0;

    all_svals = cat(dim, all_svals, svals);
  end

  return;
end
