function splines = shapes2splines(shapes)
  
  if (~isstruct(shapes))
    shapes = struct('path', shapes);
  end

  [ngroups, nframes] = size(shapes);
  for i=1:ngroups
    for j=1:nframes
      splines(i,j) = create_spline(shapes(i,j).path);
    end
  end

  return;
end
