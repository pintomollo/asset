function weight = weight_symmetry(img, params, repeat)
  
  if (nargin == 2)
    repeat = true;
  end

  img = imnorm(img);

  alpha = params.alpha;
  beta = 20 * params.beta;
  gamma = params.gamma;

  [nrows, npts] = size(img);
  corrs = symmetric_corrcoef_mex(img.').';

  corrs = (corrs + 1) / 2;
  intens = imfilter((1 - img), fspecial('gaussian', [1 ceil(2*beta)], beta), 'symmetric');
  
  weight = alpha*corrs + (1-alpha)*intens;

  [a,s] = mymean(weight, 2);
  %weight = bsxfun(@times, weight, 1-a);
  %weight = bsxfun(@times, weight, s);
  imfs = emdc([], s-a);
  factors = gamma*imnorm(imfs(end, :)) + (1-gamma);
  weight = bsxfun(@times, weight, factors.');

  weight(isnan(weight)) = Inf;

  return;
end
