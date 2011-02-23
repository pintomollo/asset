function [entropy, mean_stddev] = posterior_decoding(emission, trans, beta, gamma)

  [nindx, nhood, npts] = size(emission);
  %nindx = length(path);

  half = floor(nhood/2);

  if (beta == 0)
    emission = ones(size(emission));
  else
    emission = log(emission) * beta;
    emission = exp(emission - repmat(max(max(emission,[],3),[],2),[1 nhood npts]));
  end
  %emission = emission ./ repmat(sum(sum(emission,3),2),[1 nhood npts]);
  if (gamma == 0)
    trans = ones(size(trans));
  else
    trans = log(trans) * gamma;
    trans = exp(trans - max(trans));
  end
  %trans(2:end-1) = trans(2:end-1) / sum(trans(2:end-1));

  forward = zeros(nindx, npts); 
  backward = zeros(nindx, npts); 
  fscale = ones(nindx,1);
  bscale = ones(nindx,1);

  forward(1,:) = shiftdim(emission(1,half+1,:),1) .* trans(1,1);
  fscale(1,1) = sum(forward(1,:));

  if (fscale(1,1) ~= 0)
    forward(1,:) = forward(1,:) / fscale(1,1);
  end

  for i=2:nindx
    tmpforw = forward(i-1, [end-half+1:end 1:end 1:half]);
    for j=1:nhood
      forward(i,:) = forward(i,:) + tmpforw(j:j+npts-1).*shiftdim(emission(i,j,:),1).*trans(1,j+1);
    end
    fscale(i,1) = sum(forward(i,:));
    if (fscale(i,1) ~= 0)
      forward(i,:) = forward(i,:) / fscale(i,1);
    end
  end

  backward(end,:) = trans(1,end);
  bscale(end,1) = sum(backward(end,:));

  if (bscale(end,1) ~= 0)
    backward(end,:) = backward(end,:) / bscale(end,1);
  end
  for i=nindx-1:-1:1 
    tmpback = backward(i+1, [end-half+1:end 1:end 1:half]);
    tmpemm = shiftdim(emission(i+1,:, [end-half+1:end 1:end 1:half]),1);
    for j=1:nhood
      backward(i,:) = backward(i,:) + tmpemm(j,end-npts-j+2:end-j+1).*trans(1,end-j).*tmpback(end-npts-j+2:end-j+1);
    end
    bscale(i,1) = sum(backward(i,:));
    if (bscale(i,1) ~= 0)
      backward(i,:) = backward(i,:) / bscale(i,1);
    end
  end

  fscale = cumsum(log(fscale));
  bscale = flipud(cumsum(flipud(log(bscale))));
  
  logPx = log(sum(forward(end,:).*trans(1,end))) + fscale(end,1);
  %logPx = log(sum(backward(1,:).*trans(1,1).*shiftdim(emission(1,half+1,:),1))) + bscale(1,1)

  forward = log(forward) + repmat(fscale,1,npts);
  backward = log(backward) + repmat(bscale,1,npts);

  %%%% Works correctly but unused
  if (nargout > 1)

    Plambda = exp(forward + backward - logPx);
    dist = repmat([1:npts],nindx,1);
    means = sum(dist .* Plambda, 2);
    mean_stds = abs(sqrt(sum((dist.^2) .* Plambda, 2) - means.^2));
    
    %implot(Plambda)
    %hold on;
    %plot(means,[1:nindx],'k');
    %plot(means - mean_stds,[1:nindx],'m');
    %plot(means + mean_stds,[1:nindx],'m');

    mean_stddev = mean(mean_stds);

    entropy = Plambda;

    %figure;
    %implot(entropy)
    %figure;
    %plot(mean_stds)

    return;
  end

  emission = log(emission);
  trans = log(trans);

  entropy = 0;

  forward_spec = forward(1:end-1,:);
  forward_spec = forward_spec(:,[end-half+1:end 1:end 1:half]);

  backward = backward(2:end,:);
  emission = emission(2:end,:,:);

  % Working but unused
  %full_joint = zeros(nindx-1,nhood,npts);
  full_entr = zeros(nindx-1,nhood,npts);

  for i=1:nhood
    % Computes the full joint prob.
    %full_joint(:,i,:) = exp(forward_spec(:,i:npts+i-1) + squeeze(emission(:,i,:)) + trans(i+1) + backward - logPx);

    full_entr(:,i,:) = exp(forward_spec(:,i:npts+i-1) + squeeze(emission(:,i,:)) + trans(i+1) + backward - logPx) .* (trans(i+1) + squeeze(emission(:,i,:)));
  end

  %%% Hmmm might have to improve : 0 (joint prob.) * Inf (emmission) = NaN
  full_entr(isnan(full_entr)) = 0;

  entropy = logPx - sum(sum(sum(full_entr,3),2),1);

  return;
end
