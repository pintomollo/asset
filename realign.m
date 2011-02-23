function aligned = realign(original, new_size, center, orientation)

  rotmat = [cos(orientation) sin(orientation); -sin(orientation) cos(orientation)];
  new_center = new_size(:) / 2;
  new_center = new_center([2 1],1);

  if (isstruct(original))
    rotmat = [cos(orientation) -sin(orientation); sin(orientation) cos(orientation)];

    if (original.dim > 2)
      center = [center; zeros(original.dim - 2, 1)];
      new_center = [new_center; zeros(original.dim - 2, 1)];
      rotmat = [rotmat zeros(2); zeros(2) rotmat];
    end

    aligned = fncmb(original, '-', center);
    aligned = fncmb(aligned, rotmat);
    aligned = fncmb(aligned, '+', new_center);
  else
    [n,m] = size(original);
    if (n == 2 & m > 2)
      aligned = original - repmat(center(:), 1, m);
      aligned = (aligned.' * rotmat).';
      aligned = aligned + repmat(new_center, 1, m);
    elseif (m == 2)
      aligned = original - repmat(center(:).', n, 1);
      aligned = (aligned * rotmat);
      aligned = aligned + repmat(new_center.', n, 1);
    else
      aligned = imalign(original, new_size, center, orientation);
    end
  end

  return;
end
