function [surface, volume] = spheroidal_cap(axes_length, height)

  axes_length = sort(axes_length, 'descend');

  a = axes_length(1,:);
  b = axes_length(2,:);
  c = axes_length(3,:);

  npts = numel(a);
  if (npts > 1)
    surface = NaN(1, npts);
    volume = NaN(1, npts);

    for i=1:npts
      [surface(i), volume(i)] = spheroidal_cap(axes_length(:,i), height(i));
    end

    return;
  end

  delta = 1 - c^2 / a^2;
  epsilon = 1 - c^2 / b^2;

  surface = 4 * a * b * integral(@int_func, height/a, 1);

  height = a - height;

  volume = pi * b * c * (height^2) * (3*a - height) / (3 * (a^2));

  return;

  function int = int_func(s)

    part = 1 - delta * s.^2;
    m = epsilon * (1 - s.^2) ./ part;
    [K, E] = ellipke(m);
    int = sqrt(part) .* E;

    return;
  end
end
