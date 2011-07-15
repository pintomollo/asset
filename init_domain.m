function weight = init_domain(img, params)
  
  alpha = 1 / params.alpha;
  beta = 1 / params.beta;
  
  npts = size(img, 2);
  
  init_weight = bsxfun(@minus, [1:npts], [1:npts].')+1;
  is_inverted = (init_weight <= 0);
  init_weight(is_inverted) = init_weight(is_inverted) + npts;
  
  weight = 1 ./ (1 + exp(-alpha * (init_weight - beta)));

  return;
end
