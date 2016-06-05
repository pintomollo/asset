function [center, axes_length, orientation] = fit_ellipse(pts)
%
% To convert this vector to the geometric parameters I use 
% the standard formulas, e.g., (19) - (24) in Wolfram Mathworld:
%     http://mathworld.wolfram.com/Ellipse.html
%
  center = NaN;
  axes_length = NaN;
  orientation = NaN;

  if (isempty(pts))
    return;
  end

  ndims = size(pts, 2);
  if (ndims < 2)
    return;
  end

  center = NaN(ndims, 1);
  axes_length = NaN(ndims, 1);

  goods = all(isfinite(pts), 2);
  if (~any(goods))
    return;
  end

  pts = pts(goods, :);

  pts = unique(pts, 'rows');
  if (size(pts, 1) < 10)

    return;
  else
    repeats = bsxfun(@eq, pts, pts(1,:));

    if (any(all(repeats, 1)))
      return;
    end
  end

  if (ndims == 2)
    A = ellipsefit_direct(pts(:,1), pts(:,2));

    % Correct the scaling factor of 2 present in Wolfram
    a = A(1); b = A(2) / 2; c = A(3); d = A(4) / 2; f = A(5) / 2; g = A(6);

    delta = b^2 - a*c;
    x0 = (c*d - b*f) / delta;
    y0 = (a*f - b*d) / delta;
    center = real([x0; y0]);

    numer = 2*(a*(f^2) + c*(d^2) + g*(b^2) - 2*b*d*f - a*c*g);
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

    axes_length = [a0; b0];
  else
    [centers, axes_length, R] = ellipsoidfit_direct(pts(:,1), pts(:,2), pts(:,3));

    keyboard
  end

  % Keep only the real part of the parameters, in case something went wrong
  axes_length = real(axes_length);
  orientation = real(orientation + 2*pi * (orientation < 0) - 2*pi * (orientation > 2*pi));

  return;
end
