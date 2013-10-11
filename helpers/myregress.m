function myregress(x, y, c)

  if (nargin < 3)
    c = 'b';
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

  %% From Finn J.D., 1974
  sigma = r'*r / ndof;
  C = pinv(x'*x)*sigma^2;
  ss_all = sqrt(diag(C));
  t_all = b ./ ss_all;
  p_all = 2*(1 - tcdf(abs(t_all), ndof));
  %pred_all = sigma*sqrt(full_X*C*full_X');

  %%% Directly from Wikipedia
  ss_x = sum((pos - mean(pos)).^2);
  s_slope = sum(r.^2) / (ndof * ss_x);
  %s_inter = s_slope*sum(pos.^2) / N;
  s_inter = sum(r.^2) * (1/N + mean(pos)^2 / ss_x) / ndof;
  s_predi = s_slope.*(ss_x/N + (X - mean(pos)).^2);

  t_n = tinv(0.975, ndof);

  %{
  if (has_intercept)
    T = [b(1) / sqrt(s_inter) b(2) / sqrt(s_slope)];
  else
    T = b(1) / sqrt(s_slope);
  end

  p_val = 2*(1 - tcdf(abs(T), ndof));
  %}

  %[b - t_n*sqrt([s_inter; s_slope]) b + t_n*sqrt([s_inter; s_slope])]
  %bint

  %SS_tot = sum((y - mean(y)).^2);
  %SS_err = sum((x*b - y).^2);

  %Rsq = 1 - SS_err/SS_tot

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
      scatter(pos(i,:), y(i,:), 'MarkerEdgeColor', c(i,:));
    end
    for i=1:size(y_out,1)
      scatter(pos_out(i,:), y_out(i,:), '+',  'MarkerEdgeColor', c_out(i,:));
    end
  else
    scatter(pos, y, c);
    scatter(pos_out, y_out, ['+' c]);
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
