function weight = init_domain(img, params)
  
  alpha = params.alpha;
  beta = 1 / params.beta;
  gamma = params.gamma;
  
  npts = size(img, 2);
  nsize = floor(npts / 2);
  
  init_weight = repmat([0:nsize], npts, 1); 

  %init_weight = bsxfun(@minus, [1:npts], [1:npts].')+1;
  %is_inverted = (init_weight <= 0);
  %init_weight(is_inverted) = init_weight(is_inverted) + npts;
  
  weight = 1 ./ (1 + exp(-alpha * (init_weight - beta)));
  weight(weight > 1 - gamma) = Inf;

  return;
end
