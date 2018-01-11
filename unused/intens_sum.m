function weight = intens_sum(img, params)

  alpha = params.alpha;
  beta = params.beta;
  gamma = params.gamma;
  egg = params.path;

  npts = size(img,2);

  if (isempty(egg))
    outside = false(size(img));
  else
    pos = [1:npts];
    pos = repmat(pos,size(img,1),1) ./ repmat(egg(:),1,npts);
    outside = (pos>1);
  end
  %pos = 1./(1+exp(10*(pos - 0.5)));

  weight = img;
  weight(outside) = 0;
  weight(isnan(weight)) = 0;
  weight = cumsum(weight,2);
  weight = 1 - (weight ./ repmat(weight(:,end),1,npts));

  %figure;implot(weight-beta)
  %figure;imagesc(abs(weight-beta))
  %figure;imagesc(abs(img - gamma))

  %keyboard

  %weight = (1 - img) * alpha + (1 - alpha) * weight;
  weight = abs(img - gamma) * alpha + (1 - alpha) * abs(weight - beta);
  weight(outside) = Inf;
  weight(isnan(weight)) = Inf;
  weight(isnan(img)) = Inf;

  %figure;imagesc(weight);

  %weight = weight.^0.1;

  return;
end
