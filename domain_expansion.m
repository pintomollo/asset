function fraction = domain_expansion(domain, center, cytok)

  figure;imagesc(domain);

  [nframes, npos] = size(domain);
  domain = imfilter(domain, fspecial('average', [5, 3]), 'replicate');

  figure;imagesc(domain);

  fraction = NaN(nframes, 1);

  valids = any(isnan(domain), 1);
  tol = 3;

  first = center - find(valids(1:center), 1, 'last');
  last = find(valids(center:end), 1, 'first') - 1;
  boundary = min(first, last)-1;

  domain = domain(:,[0:boundary]+center) + domain(:,[0:-1:-boundary]+center);
  domain = imadjust(imnorm(domain(1:cytok, :)));
  %intens = cumsum(domain, 2);
  %intens = bsxfun(@rdivide, intens, intens(:, end));

  figure;imagesc(domain);

  thresh = 1.1*graythresh(domain([end-5:end], :));
  valids = (domain >= thresh);
  indxs = [1:npos];

  figure;imagesc(valids);

  init_pos = max(indxs((mean(valids([end-5:end], :)) > 0.5)));
  min_pos = init_pos;

  for i = cytok:-1:1
    tmp_valids = valids(i, :);
    tmp_valids(min_pos+tol+1:end) = false;
    new_min = max(indxs(tmp_valids));
    if (isempty(new_min))
      new_min = 0;
    end
    if (new_min < min_pos-tol)
      new_min = min_pos-tol;
    end

    if (new_min < min_pos)
      min_pos = new_min;
    end

    fraction(i) = new_min;
  end
  
  figure;imagesc(domain);
  hold on;
  plot(fraction, 1:nframes, 'k');

  fraction = fraction / init_pos;

  return;
end
