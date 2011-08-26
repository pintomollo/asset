function [bests, alpha] = estimate_gamma(x, y)

  if (numel(x) == 0)
    bests = NaN;
    alpha = NaN;

    return;
  end

  orig_y = y;
  orig_x = x;

  valids = (~isnan(x) & ~isnan(y));
  x = x(valids);
  y = y(valids);

  y = y - min(y);
  
  y = y / sum(y);

  %orig_y = y;
  %imfs = emdc([], y);
  %y = sum(imfs(2:end,:), 1).';

  center = sum(x .* y);

  s = log(sum(x .* y)) - sum(log(x) .* y);
  k = 3 - s + sqrt((s - 3)^2  + 24*s) / (12*s);
  theta = center / k;

  lambda = 1/theta;
  alpha = k;

  %x = (x - center).^2;
  %sigma = sum(y .* x);

  %alpha = center.^2 / sigma;
  %lambda = center / sigma;

  %gammas = (lambda^alpha / gamma(alpha)) * x .^ (alpha - 1) .* exp(-lambda * orig_x);
  bests = nlinfit(x, y, @gamma_function, [lambda alpha]);

  %figure;
  %plot(x, y);
  %hold on;
  %plot(x, gamma_function([lambda, alpha], x), 'r')
  %plot(x, gamma_function(bests, x), 'k')

  if (nargout == 2)
    alpha = bests(2);
    bests = bests(1);
  end

  %keyboard
  return;

  pts = pts(:);
  prob = prob(:);
  npts = numel(pts);

  prob = prob / sum(prob);

  s = log(mean(pts))
  k = 0;

  [m,s] = mymean()

  return;
end

function y = gamma_function(params, x)

  if (any(~isreal(params))|any(params <= 0))
    y = 1 + exp(x) * sum(abs(real(params)));

    return;
  end

  lambda = params(1);
  alpha = params(2);

  y = (lambda^alpha / gamma(alpha)) * (x.^(alpha - 1)) .* exp(-lambda * x);

  return;
end
