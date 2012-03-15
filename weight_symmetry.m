function weight = weight_symmetry(img, params, repeat)
  
  if (nargin == 2)
    repeat = true;
  end

  alpha = params.alpha;
  beta = 20 * params.beta;

  [nrows, npts] = size(img);
  corrs = symmetric_corrcoef_mex(img.').';

  corrs = (corrs + 1) / 2;
  intens = imfilter((1 - img), fspecial('gaussian', [1 ceil(2*beta)], beta), 'symmetric');

  weight = alpha*corrs + (1-alpha)*intens;
  weight(isnan(weight)) = Inf;

  return;
end
