function myregress(x, y, c, s)

  if (nargin < 3)
    c = 'b';
    s = 'o+';
  elseif (nargin < 4)
    s = 'o+';
  end

  if (length(s) < 2)
    s(2) = '+';
  end

  npts = 100;

  y = y(:);
  N = length(y);
  if (size(x,2) == N)
    x = x.';
  end

  has_intercept = false;
  ndims = size(x,2);

  if (all(x(:,1) == 1))
    pos = x(:,2);
    has_intercept = true;
  else
    pos = x(:,1);
  end
  ndof = size(x,1) - ndims;

  x_min = min(pos)*0.95;
  x_max = max(pos)*1.05;
  X = [x_min:(x_max-x_min)/npts:x_max].';

  [b, bint, r, rint, stats] = regress(y, x);
  goods = (rint(:,1)<0 & rint(:,2)>=0);

  y_out = y(~goods);
  x_out = x(~goods,:);
  pos_out = pos(~goods);

  y = y(goods);
  x = x(goods,:);
  pos = pos(goods);
  [b, bint, r, rint, stats] = regress(y, x);

  full_X = bsxfun(@power, X, [1:ndims]-has_intercept);
  Y = sum(bsxfun(@times, full_X, b(:).'), 2);

  ndof = size(x,1) - ndims;

  %% From Finn J.D., 1974
  %p. 100, Methods to unbias estimator of sigma
  sigma = sqrt(r'*r / ndof);
  %p. 100 covariance matrix using the variance-covariance factors
  C = inv(x'*x)*sigma^2;
  %p. 100 the standard error of the parameters
  ss_all = sqrt(diag(C));
  %p. 101 95\% confidence interval, exact same result as bint !
  %conf_int = ss_all*tinv(0.975, ndof);
  %conf_int = [-conf_int conf_int] + b(:, [1 1]);
  %p. 101, transforming into a t-distribution around 0 (b-0)
  t_all = b ./ ss_all;
  % compute the corresponding p-value of the interval
  p_all = 2*tcdf(-abs(t_all), ndof);

  %%% Directly from Wikipedia
  ss_x = sum((pos - mean(pos)).^2);
  s_slope = sum(r.^2) / (ndof * ss_x);
  %s_inter = s_slope*sum(pos.^2) / N;
  s_inter = sum(r.^2) * (1/N + mean(pos)^2 / ss_x) / ndof;
  s_predi = s_slope.*(ss_x/N + (X - mean(pos)).^2);

  t_n = tinv(0.975, ndof);

  %figure;
  hold on;
  ncolors = size(c,1);
  if (ncolors>1)
    if (ncolors ~= N)
      warning('Missing colors for several data points');
      c(end+1:N,:) = repmat(c(end,:), N-ncolors, 1);
    end
    c_out = c(~goods, :);
    c = c(goods, :);
    for i=1:size(y,1)
      scatter(pos(i,:), y(i,:), s(1), 'MarkerEdgeColor', c(i,:));
    end
    for i=1:size(y_out,1)
      scatter(pos_out(i,:), y_out(i,:), s(2),  'MarkerEdgeColor', c_out(i,:));
    end
  else
    scatter(pos, y, s(1), 'MarkerEdgeColor', c);
    scatter(pos_out, y_out, s(2), '+', 'MarkerEdgeColor', c);
  end
  plot(X, Y, 'k');
  plot(X, Y+t_n*sqrt(s_predi), 'k');
  plot(X, Y-t_n*sqrt(s_predi), 'k');

  y_limits = ylim;
  x_limits = xlim;

  range_x = diff(x_limits);
  range_y = diff(y_limits);

  text(x_limits(1) + [0.05;0.05;ones(ndims,1)*0.05]*range_x, y_limits(end) - [0.1;0.2;0.2+(([1:ndims].')/10)]*range_y, {['R^2 = ' num2str(stats(1))]; ['p_{val} = ' num2str(p_all(:).')]; num2str([b bint])});

  return;
end
