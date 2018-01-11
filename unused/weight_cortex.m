function weight = weight_cortex(img, params)

  alpha = params.alpha;
  beta = params.beta;
  gamma = 1/params.gamma;
  delta = params.delta;
  epsilon = params.epsilon;
  
  inner = params.path(:,1);
  egg = params.path(:,2);
  outer = params.path(:,3);

  [nrows,npts] = size(img);
  
  pos = [1:npts];
  pos = repmat(pos,nrows,1) ./ repmat(inner,1,npts);
  outside = (pos>1);

  edges = imadm_mex(img);
  
  gap = double(edges==0);

  edges(outside) = 0;
  %figure;imshow(edges,[],'InitialMagnification','fit');

  %entr = imnorm(entropyfilt(img, ones(1,9)));
  %entr(isnan(entr)) = 0;
  %entr = 1 - (1 - entr).^gamma;
  %figure;imshow(entr,[],'InitialMagnification','fit');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  entr = edges.^gamma;
  edges = (1-edges).^gamma;

  entr(isnan(entr)) = 0;
  entr = cumsum(entr, 2);
  entr = 1 - (entr ./ repmat(entr(:,end),1,npts));
  %figure;imshow(entr,[],'InitialMagnification','fit');

  intens = abs(img - mean(img(outside & ~isnan(img))));
  intens = (1 - imnorm(intens)).^gamma;
  %figure;imshow(intens,[],'InitialMagnification','fit');

  egg_intens = bilinear_mex(img, egg, [1:nrows].');
  egg_intens(isnan(egg_intens)) = 0;
  egg_intens = imnorm(abs(img - repmat(egg_intens,1,npts)));

  weight = (edges * delta + gap * (1 - delta)) * alpha + (entr * beta + (intens * epsilon + egg_intens * (1 - epsilon)) * (1 - beta))* (1 - alpha);
  %weight = (edges * delta + gap * (1 - delta)) * alpha + (entr * beta + intens * (1 - beta))* (1 - alpha);
  weight(isnan(weight)) = Inf;
  weight(isnan(img)) = Inf;
  weight(outside) = Inf;
  %figure;imshow(weight,[],'InitialMagnification','fit');

  return;
end
