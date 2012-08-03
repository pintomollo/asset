function weight = weight_expansion(img, params)

  alpha = params.alpha;
  beta = 2*params.beta;

  thresh = beta*graythresh(img([end-5:end], :));
  intens = imnorm(abs(imnorm(img) - thresh));

  pos = cumsum(1-imnorm(img), 2);
  pos = bsxfun(@rdivide, pos, pos(:, end));

  weight = alpha*intens + (1-alpha)*pos;

  return;
end
