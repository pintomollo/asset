function [H, pval] = myttest(values, indexes, tails)

  if (nargin == 1)
    indexes = repmat([1:size(values,2)], size(values,1), 1);
    indexes = indexes(:);
    values = values(:);
    tails = 'both';
  elseif (nargin == 2)
    tails = 'both';
  end

  groups = unique(indexes(:)).';

  if (numel(groups) == 1)
    [H,pval] = ttest(values(:));

    return;
  end

  ngroups = length(groups);
  H = zeros(ngroups);
  pval = H;

  for i = 1:ngroups
    for j = 1:ngroups
      [H(i,j), pval(i,j)] = ttest2(values(indexes==groups(i)), values(indexes==groups(j)), 0.05, tails);
    end
  end

  return;
end
