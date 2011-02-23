function weight = intens_filt(img, params)

  alpha = params.alpha;

  weight = imfilter(img, params.filt, 'symmetric');
  weight = imnorm(weight);
  weight(isnan(weight)) = Inf;

  %figure;imshow(weight,[],'InitialMagnification','fit');

  weight = img * alpha + (1 - alpha) * weight;

  %figure;imshow(weight,[],'InitialMagnification','fit');

  return;
end
