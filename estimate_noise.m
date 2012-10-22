function [gamma] = estimate_noise(img)

  mfilter = 3;
  block_size = 21;
  middle = ((block_size-1)/2)+1;
  bin_size = 10;

  row = -ones(1, block_size);

  filter = diag(row);
  filter = filter + filter.';
  filter(middle, :) = 3*row;
  filter(:, middle) = 3*row;
  filter(middle, middle) = 8*(block_size - 1);

  filter = filter(:);
  nelems = numel(filter);

%% Filters added ultimately
%  hfilter = filter;
%  hfilter(middle, :) = row;
%  vfilter = hfilter.';
%  pdfilter = diag(row);
%  ndfilter = pdfilter.';
%  ldfilter = filter;
%  ldfilter

  noisefree = median_mex(img, mfilter);
  noisy = img - noisefree;

  blocks = blockproc(img, [block_size block_size], @analyze_block);
  blocks = reshape(blocks, [], 3);

  blocks = sortrows(blocks);
  target_var = mean(blocks(1:3, 3));

  goods = (blocks(:,3) < 3*target_var);
  blocks = blocks(goods, :);

  gauss_noise = mean(blocks(1:20, 2:3));

  noisefree = noisefree(:);
  noisy = noisy(:);
  edges = [min(noisefree)-1:bin_size:max(noisefree)+1];
  [counts, map] = histc(noisefree, edges);
  nbins = length(counts);

  %%%%%%%% Use abacus stuff instead...
  minsize = 84;

  vars = NaN(nbins, 1);
  for i=1:nbins
    if (counts(i) > minsize)
      vars(i) = std(noisy(map==i))^2;
    end
  end

  goods = isfinite(vars);
  vars = vars(goods);
  edges = edges([goods; false]);

  keyboard

  return;

  function res = analyze_block(block_struct)
    
    data = block_struct.data(:);

    if (numel(data) == nelems)
      homog = sum(data.*filter);
      [means, stds] = mymean(data);

      res = cat(3, cat(3, homog, means), stds^2);
    else
      res = cat(3, cat(3, Inf, 0), 0);
    end

    return;
  end
end
