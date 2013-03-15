function weight = weight_expansion(img, params)

  if (size(img, 1) <= 10)
    warning('Not enough data points to estimate the weight matrix properly');
    weight = ones(size(img));
    
    return;
  end

  minmax = prctile(img(isfinite(img)), [0.1 99.9]);
  img = imnorm(img, minmax(1), minmax(2));
  %img = imnorm(img);

  bads = isnan(img);
  img(bads) = 0;

  alpha = params.alpha;
  beta = 2*params.beta;
  gamma = params.gamma;

  pos = cumsum(img, 2);
  pos = imnorm(bsxfun(@minus, pos(:, end), pos));
  thresh = graythresh(pos([end-5:end], :));
  inners = any(pos(end-5:end, :) > thresh, 1);
  inners(end) = false;

  thresh = graythresh(img([end-5:end], :));
  intens = imnorm(abs(img - beta*thresh));

  %pos = cumsum(1-imnorm(img), 2);
  %pos = bsxfun(@rdivide, pos, pos(:, end));

  slopex = differentiator(img, 2, 19);
  slopey = differentiator(img, 1, 19);
  slope = imnorm(min(slopex, -slopey));

  weight = alpha*(gamma*intens + (1-gamma)*pos) + (1-alpha)*slope;
  weight(bads) = Inf;
  weight(end-5:end, inners) = Inf;

  %figure;imagesc(weight)
  %keyboard

  return;
end
