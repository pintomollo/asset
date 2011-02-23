function spline = elliptic2spline(spline, center, axes_length, orientation)

  rotmat = [[cos(orientation) -sin(orientation)] .* axes_length'; ...
            [sin(orientation) cos(orientation)] .* axes_length'];

  spline = fncmb(fncmb(spline,rotmat),center);
  
  return;
end
