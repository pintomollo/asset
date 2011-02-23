function spline = elliptic_spline(pts, center, axes_length, orientation, sort_pts)

  if (nargin < 5)
    sort_pts = false;
  end

  ell_pts = carth2elliptic(pts, center, axes_length, orientation);

  if (sort_pts)
    [junk, indxs] = sort(ell_pts(:,1));

    ell_pts = ell_pts(indxs,:);
    ell_pts = ell_pts([end-2:end 1:end 1:3],:);
    ell_pts(1:3,1) = ell_pts(1:3,1) - 2 * pi;
    ell_pts(end-2:end,1) = ell_pts(end-2:end,1) + 2 * pi;
  end

  spline = create_spline(ell_pts(:,2),ell_pts(:,1));

  return;
end
