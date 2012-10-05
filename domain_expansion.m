function fraction = domain_expansion(domain, center, cytok, opts, domain_half)

  if (nargin == 5)
    [domain, domain_half, center, cytok, opts] = deal(domain, center, cytok, opts, domain_half);

    goods = ~isnan(domain);
    goods_half = ~isnan(domain_half);

    domain(~goods) = 0;
    domain_half(~goods_half) = 0;

    domain = (domain + domain_half) ./ (goods + goods_half);
  end

%  figure;imagesc(domain);
  bads = isnan(domain);
  domain(bads) = 0;


  [nframes, npos] = size(domain);
  domain = imfilter(domain, fspecial('average', [5, 3]), 'replicate');

%  figure;imagesc(domain);

  fraction = NaN(nframes, 1);

  valids = any(isnan(domain), 1);
  tol = 3;

  first = center - find(valids(1:center), 1, 'last');
  if (isempty(first))
    first = center-1;
  end
  last = find(valids(center:end), 1, 'first') - 1;
  if (isempty(last))
    last = length(valids)-center;
  end
  boundary = min(first, last)-1;

  if (boundary <= 10)
    if (max(first, last) > 10)
      if (first > last)
        domain = domain(:,[0:-1:-first]+center);
      else
        domain = domain(:,[0:last]+center);
      end
    else
      warning('Detected width of the domain is too small for analysis');

      return;
    end
  else
    domain = domain(:,[0:boundary]+center) + domain(:,[0:-1:-boundary]+center);
  end

  domain = imnorm(domain(1:cytok, :));
  path = dynamic_programming(domain, opts.segmentation_parameters.domain_expansion.cortex_params, opts.segmentation_parameters.domain_expansion.scoring_func, opts.segmentation_parameters.domain_expansion.cortex_weights, opts);

  fraction = path / max(path(end-5:end));

  %{
  domain = imadjust(imnorm(domain(1:cytok, :)));

  params = get_struct('smoothness_parameters');
  weights = get_struct('data_parameters');
  opts = get_struct('ASSET');
  opts.force_circularity = false;
  opts.dp_method = 'normal';

  weights.alpha = 0.55;
  weights.beta = 0.35;
  weights.gamma = 0.75;

  params.init = 1;
  params.nhood = 21;
  params.alpha = 0.65;
  params.beta = 0.75;
  params.gamma = 0.15;
  params.prohibit = 'horiz';
  params.spawn_percentile = 0.1;

  path = dynamic_programming(domain, params, @weight_expansion, weights, opts);
  %}

  %figure;imagesc(domain);
  %axis([1 size(domain, 2) 1 size(domain, 1)]);
  %hold on;plot(path, [1:size(domain, 1)], 'k')

  %keyboard

  %{
  thresh = 1.1*graythresh(domain([end-5:end], :));
  valids = (domain >= thresh);
  indxs = [1:npos];

%  figure;imagesc(valids);

  init_pos = max(indxs((mean(valids([end-5:end], :)) > 0.5)));
  
  if (isempty(init_pos))
    return;
  end
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
  
%  figure;imagesc(domain);
%  hold on;
%  plot(fraction, 1:nframes, 'k');

  fraction = fraction / init_pos;
  %}

  return;
end
