function [fraction, max_width, cell_width, raw_domain, pos, path_center] = domain_expansion(domain, center, cytok, opts, domain_half)

  mymovie = [];
  raw_domain = [];
  pos = [];
  if (nargin < 2)
    error('Minimal number of arguments is 2');
  elseif (nargin == 2)
    mymovie = domain;
    opts = center;

    [raw_domain, ruffles, theta] = gather_quantification(mymovie, opts);

  bads = isnan(raw_domain);
  domain = inpaint_nans(raw_domain);
  domain = gaussian_mex(domain, 0.67);
  domain(bads) = NaN;

    domain = imnorm(domain);
    opts = load_parameters(opts, 'domain_center.txt');
    opts.quantification.weights.filt = ruffles;
    opts.quantification.params.init = (1-exp(-theta.^2/(2*(opts.quantification.params.spawn_percentile(1)/10)^2)));
    params = opts.quantification;

    path_center = dynamic_programming(domain, opts.quantification.params, opts.quantification.scoring_func, opts.quantification.weights, opts);

    if (opts.verbosity == 3)

      center = find(theta == 0);
      pos_indx = [center:50:size(domain,2)];
      tmp = [center:-50:1];
      pos_indx = [tmp(end:-1:2) pos_indx];

      figure;imagesc(ruffles);

      figure;subplot(2,1,1)
      imagesc(domain);
      set(gca, 'XTick', pos_indx, 'XTickLabel', theta(pos_indx));
      hold on;
      plot(path_center, 1:length(path_center));
      subplot(2,1,2);
      plot(opts.quantification.params.init);
      xlim([1 size(domain, 2)]);

      figure;imagesc(opts.quantification.scoring_func(domain, opts.quantification.weights));
    end

    [domain, ruffles, pos, indx] = align_domain(domain, ruffles, path_center, opts);
    [domain, boundary] = crop_domain(domain, indx);
    pos = pos([-boundary:boundary]+indx);

    %raw_domain = domain;

    opts = load_parameters(opts, 'domain_expansion.txt');
    time = get_manual_timing(mymovie, opts);
    cytok = time(end);
    center = boundary + 1;

    raw_domain = domain(1:cytok, :);
    domain = nanmean(cat(3, domain(:, center:end), domain(:,center:-1:1)), 3);

  elseif (nargin == 5)
    [domain, domain_half, center, cytok, opts] = deal(domain, center, cytok, opts, domain_half);

%    goods = ~isnan(domain);
%    goods_half = ~isnan(domain_half);

%    domain(~goods) = 0;
%    domain_half(~goods_half) = 0;

%    raw_domain = (domain + domain_half) ./ (goods + goods_half);

    domain = nanmean(cat(3, domain, domain_half), 3);

    if (center == size(domain,2))
      domain = domain(:,end:-1:1);
    elseif (center ~= 1)
      domain = nanmean(cat(3, domain(:, center:end), domain(:,center:-1:1)), 3);
    end
    boundary = size(domain, 2);
  else
    if (center == size(domain,2))
      domain = domain(:,end:-1:1);
    elseif (center ~= 1)
      domain = nanmean(cat(3, domain(:, center:end), domain(:,center:-1:1)), 3);
    end
    boundary = size(domain, 2);
  end

%  figure;imagesc(domain);
%  domain(bads) = 0;


  [nframes, npos] = size(domain);
%  domain = imfilter(domain, fspecial('average', [5, 3]), 'replicate');

%  figure;imagesc(domain);

  fraction = NaN(nframes, 1);
  max_width = NaN;
  cell_width = NaN;

  %bads = isnan(domain);
  %domain(bads) = 0;
  %domain = imfilter(domain, fspecial('average', [5, 3]), 'replicate');
  %domain(bads) = NaN;
  %valids = any(isnan(domain(~bads,:)), 1);
  %tol = 3;

  %{
  first = center - find(valids(1:center), 1, 'last');
  if (isempty(first))
    first = center-1;
  end
  last = find(valids(center:end), 1, 'first') - 1;
  if (isempty(last))
    last = length(valids)-center;
  end
  boundary = min(first, last)-1;

  if (boundary <= 5)
    if (max(first, last) > 5)
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
  %}
%  domain = inpaint_nans(raw_domain);
%  domain = nanmean(cat(3, domain(:, center:end), domain(:,center:-1:1)), 3);
  %domain = domain(:, center:end) + domain(:,center:-1:1);
%  domain = gaussian_mex(domain, 0.67);
  %domain = domain(:,center:end) + domain(:,center:-1:1);

  %domain = gaussian_mex(domain, 1.5);
  domain = imnorm(domain(1:cytok, :));
  path = dynamic_programming(domain, opts.segmentation_parameters.domain_expansion.cortex_params, opts.segmentation_parameters.domain_expansion.scoring_func, opts.segmentation_parameters.domain_expansion.cortex_weights, opts);

    if (opts.verbosity == 3)
      pos_indx = [0:50:boundary];

      if (isempty(raw_domain))
        raw_domain = [domain(:,end:-1:2) domain];
      end

      figure;
      imagesc(raw_domain);
      if (~isempty(pos))
        set(gca, 'XTick', [-pos_indx([end:-1:2]) pos_indx]+center, 'XTickLabel', pos([-pos_indx([end:-1:2]) pos_indx]+center));
      end

      figure;imagesc(domain);
      hold on;plot(path, 1:length(path));
      if (~isempty(pos))
        set(gca, 'XTick', pos_indx+1, 'XTickLabel', pos(pos_indx+center));
      end

      figure;imagesc(opts.segmentation_parameters.domain_expansion.scoring_func(domain, opts.segmentation_parameters.domain_expansion.cortex_weights));
    end
  %figure;
  %imagesc(domain);hold on;
  %plot(path, 1:length(path), 'k');figure;

  cell_width = median(2*sum(~isnan(domain(max(end-5, 1):end, :)), 2) - 1) * opts.quantification.resolution;
  max_width = max(path(max(end-5, 1):end));
  fraction = path / max_width;

  if (opts.quantification.resolution > 0)
    max_width = max_width * opts.quantification.resolution;
  end

%  if (~isempty(mymovie))
%    mymovie.quantification.
%  end

  return;
end
