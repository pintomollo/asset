function [H, pval] = myttest(values, indexes)

  groups = unique(indexes(:)).';

  if (numel(groups) == 1)
    [H,pval] = ttest(values(:));

    return;
  end

  H = zeros(length(groups));
  pval = H;

  for i = groups
    for j = groups
      [H(i,j), pval(i,j)] = ttest2(values(indexes==i), values(indexes==j), [], 'right');
    end
  end

  return;
end
