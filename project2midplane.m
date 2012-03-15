function [midx, midy] = project2midplane(ptsx, ptsy, center, axes_length, orient, zpos, thresh)

  midx = [];
  midy = [];

  if (isempty(ptsx))
    return;
  end

  if (nargin==5)
    [center, axes_length, orient, zpos] = deal(ptsy, center, axes_length, orient);
    ptsy = [];
    thresh = 0.1;
  elseif (narign == 6)
    if (numel(ptsy) == numel(ptsx))
      thresh = 0.1;
    else
      [center, axes_length, orient, zpos, thresh] = deal(ptsy, center, axes_length, orient, zpos);
      ptsy = [];
    end
  end

  if (size(ptsx,2) > 4)
    ptsx = ptsx.';
    ptsy = ptsy.';
  end

  if (size(ptsx, 2) == 4)
    [pts_x, pts_y] = project2midplane(ptsx(:,1), ptsx(:,2), center, axes_length, orient, zpos, thresh);
    [stds_x, stds_y] = project2midplane(ptsx(:,3), ptsx(:,4), [0; 0], axes_length, orient, zpos, thresh);

    if (nargout == 2)
      midx = [pts_x stds_x];
      midy = [pts_y stds_y];
    else
      midx = [pts_x pts_y stds_x stds_y];
    end

    return;
  elseif (size(ptsx, 2) == 2)
    ptsy = ptsx(:,2);
    ptsx = ptsx(:,1);
  end

  corient = cos(orient);
  sorient = sin(orient);

  % Centering
  x = ptsx - center(1);
  y = -(ptsy - center(2));

  % Rotating
  tmp_x = x*corient + y*sorient;
  tmp_y = -x*sorient + y*corient;
  x = tmp_x;
  y = tmp_y;

  % Variables for the projection
  n = abs(axes_length(1) ./ x);
  b = zpos.^2 / axes_length(3)^2;
  m = 1 ./ sqrt(1 - b);
  Cinf = m^2;

  % Coefficients for the projection
  coef = 1 ./ (1 - b.* (n.^2 ./ (n.^2 - 1)));

  % Variables for the sigmoid
  obj = 1 / (Cinf * (1 + thresh));
  xt = sqrt((1 - obj) / (1 - b - obj));
  th = (xt + m) / 2;
  phi = (log(1 - thresh) - log(thresh)) / (xt - th);

  Ay = (Cinf*(1+thresh) - m) / (1 - thresh);
  Ax = (1 - m) / (1 - thresh);
  
  % Threshold the distance to replace
  near = (n < xt);

  coef(near) = Ay ./ (1 + exp(-phi*(n(near) - th))) + m;
  coefx = ones(size(x));
  coefx(near) = Ax ./ (1 + exp(-phi*(n(near) - th))) + m;

  y = y .* coef;
  x = x .* coefx;

  % Rotating & centering
  midx = x * corient - y*sorient + center(1);
  midy = -(x * sorient + y * corient) + center(2);

  if(nargout==1)
    midx = [midx midy];
    midy = [];
  end

  return;
end
