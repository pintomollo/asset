function [midx, midy] = project2midplane(ptsx, ptsy, center, axes_length, orient, zpos)

  midx = [];
  midy = [];

  if (isempty(ptsx))
    return;
  end

  if (nargin==5)
    [center, axes_length, orient, zpos] = deal(ptsy, center, axes_length, orient);
  end

  if (size(ptsx,2) > 4)
    ptsx = ptsx.';
  end

  if (size(ptsx, 2) == 4)
    [pts_x, pts_y] = project2midplane(ptsx(:,1), ptsx(:,2), center, axes_length, orient, zpos);
    [stds_x, stds_y] = project2midplane(ptsx(:,3), ptsx(:,4), [0; 0], axes_length, orient, zpos);

    if (nargout == 2)
      midx = [pts_x stds_x];
      midy = [pts_y stds_y];
    else
      midx = [pts_x pts_y stds_x stds_y];
    end

    return;
  end

  ptsy = ptsx(:,2);
  ptsx = ptsx(:,1);

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

  % Rescaling
  y = real(y ./ sqrt(1 - (axes_length(1)^2 * zpos^2) ./ (axes_length(2)^2 .* (axes_length(1)^2 - x.^2))));

  % Rotating & centering
  midx = x * corient - y*sorient + center(1);
  midy = -(x * sorient + y * corient) + center(2);

  if(nargout==1)
    midx = [midx midy];
    midy = [];
  end

  return;
end
