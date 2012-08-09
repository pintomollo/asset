function weight = weight_expansion(img, params)

  alpha = params.alpha;
  beta = 2*params.beta;
  gamma = params.gamma;

  thresh = beta*graythresh(img([end-5:end], :));
  intens = imnorm(abs(imnorm(img) - thresh));

  pos = cumsum(1-imnorm(img), 2);
  pos = bsxfun(@rdivide, pos, pos(:, end));

  slopex = differentiator(img, 2);
  slopey = differentiator(img, 1);
  slope = imnorm(min(slopex, -slopey));

  weight = alpha*(gamma*intens + (1-gamma)*pos) + (1-alpha)*slope;
  %figure;imagesc(weight)

  %keyboard

  return;
end
