function [y1, y2] = front(p, x)

  fliped = false;
  if (size(x,1) ~= 1)
    x = x.';
    fliped = true;
  end

  if (numel(p) == 6)
    p = [p(1:2) NaN p(3:end)];
  end

  y1 = piecewise_background(p(:, 1:end-1), x);
  y2 = gaussian([p(:, [1:2 end]) zeros(size(p(:,1)))], x);

  if (fliped)
    y1 = y1.';
    y2 = y2.';
  end

  if (nargout ~= 2)
    y1 = y1 + y2;
  end

  return;
end

function y = piecewise_background(p, x)
  
  if (size(p, 2) == 5)
    p = [p(:, [1:2]) NaN(size(p,1), 1) p(:, 3:end)];
  end

  problems = (any(p(:, 2:3) < 0, 2) | any(~isfinite(p(:,[1:2 4:end])), 2));
  invags = (isnan(p(:, 3)) & ~problems);

  y = zeros(size(p, 1), length(x));

  p(:, 2) = 1.75*p(:,2);

  neg = bsxfun(@lt, x, p(:, 1) - p(:, 2));
  pos = bsxfun(@gt, x, p(:, 1) + p(:, 2));
  middle = (~neg & ~pos);

  knots = [-p(:, 5).*(-p(:, 2)) + p(:, 6), ...
           p(:, 6) - p(:, 4).*(1 - exp(-p(:, 3).*p(:, 2)))];

  tmp = bsxfun(@minus, p(:, 6), bsxfun(@times, p(:, 4), (1 - exp(bsxfun(@times, -p(:, 3), bsxfun(@minus, x, p(:, 1)))))));

  if (any(invags))
    knots(invags, 2) = -p(invags, 4).*p(invags,2) + p(invags,6);
    tmp(invags, :) = bsxfun(@plus, bsxfun(@times, -p(invags, 4), bsxfun(@minus, x, p(invags, 1))), p(invags, 6));
  end

  y(pos) = tmp(pos);

  tmp = bsxfun(@plus, bsxfun(@times, bsxfun(@rdivide, diff(knots, [], 2), 2*p(:, 2)), bsxfun(@minus, x, p(:, 1))), mean(knots, 2));
  y(middle) = tmp(middle);

  tmp = bsxfun(@plus, bsxfun(@times, -p(:, 5), bsxfun(@minus, x, p(:, 1))), p(:, 6));
  y(neg) = tmp(neg);

  y(problems, :) = 0;

  return;
end

function y = gaussian(p, x)

  if (any(isnan(p)) | any(isinf(p)))
    y = zeros(size(x));

    return;

  elseif (size(p, 2) < 4)
    p = [p zeros(size(p, 1), 1)];
  end

  y = bsxfun(@plus, bsxfun(@times, p(:, 3), exp(bsxfun(@rdivide, -(bsxfun(@minus,x,p(:,1)).^2), (2*(p(:, 2).^2))))), p(:, 4));

  return;
end


