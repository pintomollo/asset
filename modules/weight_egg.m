function weight = weight_egg(img, params)

  alpha = params.alpha;
  beta = 1 / params.beta;

  [nrows,npts] = size(img);

  edge = imadm_mex(img);

  outside = edge.^beta;
  edge = (1-edge).^beta;

  outside(isnan(outside)) = 0;
  outside = cumsum(outside,2);
  outside = 1 - (outside ./ repmat(outside(:,end),1,npts));

  weight = edge * alpha + outside * (1 - alpha);
  weight(isnan(weight)) = Inf;
  weight(isnan(img)) = Inf;

  %figure;imshow(weight,[],'InitialMagnification','fit');

  return;
end
