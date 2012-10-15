function weight = weight_symmetry(img, params)
  
  img = imnorm(img);

  alpha = params.alpha;
  beta = 20 * params.beta;
  if (beta < 1)
    beta = 1;
  end
  gamma = params.gamma;
  delta = params.delta;

  depth = params.filt;

  if (isempty(depth))
    error('Missing invagination data for domain centration');
  end

  ruffles = depth;
  outside = ~isfinite(ruffles);
  ruffles(outside) = 0;
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

  %figure;imagesc(weight)

  return;
end
