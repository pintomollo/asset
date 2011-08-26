function [center, axes_length, orientation] = fit_ellipse(X,Y)
%
% To convert this vector to the geometric parameters I use 
% the standard formulas, e.g., (19) - (24) in Wolfram Mathworld:
%     http://mathworld.wolfram.com/Ellipse.html
%
  center = NaN(2,1);
  axes_length = center;
  orientation = NaN;

  if (isempty(X))
    return;
  end

  if (nargin == 1)
    if (size(X, 1)==2 & size(X, 2)>2)
      X = X.';
    end

    Y = X(:,2);
    X = X(:,1);
  else
    X = X(:);
    Y = Y(:);
  end

  indx = (any(isnan([X Y]),2) | any(isinf([X Y]),2));

  if (any(indx))
    if (all(indx))

      return;
    end

    X = X(~indx);
    Y = Y(~indx);
  end

  XY = unique([X Y], 'rows');
  if (size(XY, 1) < 4)

    return;
  end

  A = EllipseDirectFit(XY);
  if (isempty(A))
    return;
  end

  % Correct the scaling factor of 2 present in Wolfram
  a = A(1); b = A(2) / 2; c = A(3); d = A(4) / 2; f = A(5) / 2; g = A(6);

  delta = b^2 - a*c;
  x0 = (c*d - b*f) / delta;
  y0 = (a*f - b*d) / delta;
  center = real([x0; y0]);

  numer = 2*(a*f^2 + c*d^2 + g*b^2 - 2*b*d*f - a*c*g);
  denom = sqrt((a-c)^2 + 4*b^2);
  a0 = sqrt(numer / (delta * (denom - (a + c))));
  b0 = sqrt(numer / (delta * (-denom - (a + c))));

  if (b == 0)
    if (a < c)
      orientation = 0;
    else
      orientation = pi / 2;
    end
  else
    if (a < c)
      orientation = acot((a - c) / (2*b)) / 2;
    else
      orientation = acot((a - c) / (2*b)) / 2 + pi/2;
    end
  end
  % Correct for the inverted Y axis in images
  orientation = -orientation;

  % Reorder the axes length
  if (a0 < b0)
    tmp = b0;
    b0 = a0;
    a0 = tmp;

    orientation = orientation + pi/2;
  end

  % Keep only the real part of the parameters, in case something went wrong
  axes_length = real([a0; b0]);
  orientation = real(orientation + 2*pi * (orientation < 0) - 2*pi * (orientation > 2*pi));

  return;
end
