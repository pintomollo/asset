function [carth_pts, carth_y] = elliptic2carth(ptso, ptsr, center, axes_length, orient)
  
  if (numel(ptso)==0)
    carth_pts = ptso;

    if(nargout==2)
      carth_y = ptso;
    end

    return;
  end

  if (nargin < 5)
    if (nargin == 1)
      [center, axes_length, orient] = fit_ellipse(ptso);
    else
      [center, axes_length, orient] = deal(ptsr, center, axes_length);
    end

    if (size(ptso,2) > 4)
      ptso = ptso.';
    end

    if (size(ptso, 2) == 4)
      [pts_x, pts_y] = elliptic2carth(ptso(:,1), ptso(:,2), center, axes_length, orient);
      [stds_x, stds_y] = elliptic2carth(ptso(:,3), ptso(:,4), [0; 0], axes_length, orient);

      if (nargout == 2)
        carth_pts = [pts_x stds_x];
        carth_y = [pts_y stds_y];
      else
        carth_pts = [pts_x pts_y stds_x stds_y];
      end

      return;
    end

    ptsr = ptso(:,2);
    ptso = ptso(:,1);
  else
    ptso = ptso(:);
    ptsr = ptsr(:);
  end

  orient = orient; 
  corient = cos(orient);
  sorient = sin(orient);

  O = ptso;
  r = ptsr;

  carth_pts = zeros(length(ptso),2);
  carth_pts(:,1) = axes_length(1)*r.*cos(O)*corient - axes_length(2)*r.*sin(O)*sorient + center(1);
  carth_pts(:,2) = -(axes_length(1)*r.*cos(O)*sorient + axes_length(2)*r.*sin(O)*corient) + center(2);

  if(nargout==2)
    carth_y = carth_pts(:,2);
    carth_pts(:,2) = [];
  end

  return;
end
