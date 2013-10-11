function h = find_cap_size(axes_length, arc_length)

  a = axes_length(1,:);
  b = axes_length(2,:);
  c = axes_length(3,:);

  npts = numel(a);
  if (npts > 1)
    h = NaN(1, npts);

    for i=1:npts
      h(i) = find_cap_size(axes_length(:,i), arc_length(i));
    end

    return;
  end

  e = sqrt(1 - b.^2 / a.^2);
  h = fzero(@(u)(arc_length - get_length(u)), [-a a]);

  return;

  function curr_arc = get_length(x)

    y = b*sqrt(1 - x.^2/a.^2);
    curr_angle = atan2(y, x);
    if (curr_angle < 0)
      curr_angle = curr_angle + 2*pi;
    end
    [F, E] = elliptic12(curr_angle, e);
    curr_arc = 2 * a * E;

    return;
  end
end
