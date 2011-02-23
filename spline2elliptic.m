function spline = spline2elliptic(spline, center, axes_length, orientation)

  rotmat = [[cos(orientation) sin(orientation)] / axes_length(1); ...
            [-sin(orientation) cos(orientation)] / axes_length(2)];

  spline = fncmb(fncmb(spline,-center),rotmat);
  
  return;
end
