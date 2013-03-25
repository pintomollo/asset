function weight = weight_symmetry(img, params)
  
  minmax = prctile(img(isfinite(img)), [0.1 99.9]);
  img = imnorm(img, minmax(1), minmax(2));
  %img = imnorm(img);

  alpha = params.alpha;
  beta = 20 * params.beta;
  if (beta < 1)
    beta = 1;
  end
  gamma = 200*params.gamma;
  delta = params.delta;

  depth = params.filt;
  %center = params.path;

  if (isempty(depth))
    error('Missing invagination data for domain centration');
  end
  %if (isempty(center))
  %  error('Missing cortical position data for domain centration');
  %end

  lengths = repmat(nanstd(img, [], 2), 1, size(img,2));

  outside = (~isfinite(depth) | ~isfinite(img));

  depth(outside) = 0;
  img(outside) = 0;

  %ruffles = imfilter(depth, fspecial('gaussian', [1 ceil(4*beta)], 2*beta), 'symmetric');
  %ruffles = imfilter(ruffles, fspecial('gaussian', [ceil(beta) 1], 0.5*beta), 'symmetric');
  ruffles = median_mex(depth, ceil(beta));
  ruffles = imfilter(ruffles, fspecial('gaussian', [1 ceil(4*beta)], 2*beta), 'symmetric');

  %domain = imfilter(img, fspecial('gaussian', [1 ceil(4*beta)], 2*beta), 'symmetric');
  %domain = imfilter(domain, fspecial('gaussian', [ceil(beta) 1], 0.5*beta), 'symmetric');
  domain = median_mex(img, ceil(beta));
  domain = imfilter(domain, fspecial('gaussian', [1 ceil(4*beta)], 2*beta), 'symmetric');

  %lengths = gamma*size(img,2);
  lengths = gamma;

  %lengths = (imnorm(lengths)*0.5 + 0.5)*gamma*size(img,2);
  %lengths = ((1-imnorm(lengths))*0.25 + 0.75)*gamma*size(img,2);

  corrs = symmetric_corrcoef_mex(domain.', lengths.').';
  invags = symmetric_corrcoef_mex(ruffles.', lengths.').';

  corrs = (corrs + 1) / 2;
  invags = (invags + 1) / 2;

  %inners = abs(center) < 10;
  %bias_x = abs(center).^(1/2);
  %bias_x = bias_x / max(bias_x);
  %bias_y = [0:size(domain,1)-1].';
  %bias_y = bias_y / max(bias_y);

  %bias = bias_y * bias_x;

  %ruffles = imnorm(ruffles, [], 0.1);

  %vals = prctile(ruffles(:), [40 90]);
  %ruffles = imnorm(ruffles, 0, vals(2));
  ruffles = imnorm(ruffles, 0, 0.075);

  weight = alpha*(invags*delta + (1-delta)*ruffles) + (1-alpha)*corrs;
  %weight = alpha*(delta*invags + (1-delta)*ruffles) + (1-alpha)*bias;
  %weight = alpha*(gamma*invags + (1-gamma)*ruffles) + (1-alpha)*(delta*corrs + (1-delta)*bias);
  weight = bsxfun(@times, weight, imnorm(1-nanmean(domain, 2))/2 + 0.5);
  weight(isnan(weight) | outside) = Inf;
  %weight(end, ~inners) = Inf;

  return;

  ruffles = imfilter(ruffles, fspecial('gaussian', [1 ceil(4*beta)], 2*beta), 'symmetric');
  druffles = differentiator(ruffles, 2);
  ddruffles = differentiator(druffles, 2);

  maxs = (abs(druffles) < 0.05) & (ddruffles < -5e-5);
  invs = double(~maxs & ~outside);
  dist = NaN(size(maxs));

  for i=1:size(ruffles, 1)
    edges = ~outside(i,:);
    shift = find(edges, 1, 'first');
    
    if (isempty(shift))
      dist(i, :) = 1;
    else
      [pos, len, val] = boolean_domains(maxs(i, edges), true);
      pos = pos(~val) + shift;
      pos(1) = 1;
      len = len(~val);

      for j = 1:length(pos)-1
        dist(i, pos(j):pos(j+1)-1) = len(j);
      end
      dist(i, pos(end):end) = len(end);
    end
  end
  dist = round(dist * 0.625);

  min_dist = ceil(size(img, 2) / 40);
  dist(dist < min_dist) = min_dist;

  bads = isnan(img);
  img(bads) = 0;

  corrs = symmetric_corrcoef_mex(img.', dist.').';
  invags = symmetric_corrcoef_mex(depth.', dist.').';

  corrs = (corrs + 1) / 2;
  invags = (invags + 1) / 2;

  intens = imfilter((1 - img), fspecial('gaussian', [1 ceil(2*beta)], beta), 'symmetric');
  
  weight = alpha*(gamma*corrs + (1-gamma)*invags) + (1-alpha)*(delta*intens + (1-delta)*ruffles);

  %[a,s] = mymean(weight, 2);
  %imfs = emdc([], s-a);
  %factors = gamma*imnorm(imfs(end, :)) + (1-gamma);
  %weight = bsxfun(@times, weight, factors.');

  weight(isnan(weight) | bads) = Inf;
  %weight(all(isinf(weight), 2), :) = mymean(weight(:));

  %figure;imagesc(weight)

  return;
end
