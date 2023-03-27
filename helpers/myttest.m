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

  vfinites = isfinite(values);
  for i = 1:ngroups
    validsi = vfinites & indexes==groups(i);
    if (~any(validsi))
      H(i,:) = NaN;
      pval(i,:) = NaN;
    else
      for j = 1:ngroups
        validsj = vfinites & indexes==groups(j);
        if (~any(validsj))
          H(i,j) = NaN;
          pval(i,j) = NaN;
        else
          [H(i,j), pval(i,j)] = ttest2(values(validsi), values(validsj), 'alpha', 0.05, 'tail', tails);
        end
      end
    end
  end
  H = H + (pval < 0.01) + (pval < 0.001);

  return;
end
