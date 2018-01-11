function conc = compute_concavity(pts, thresh)

  npts = size(pts, 1);
  conc = false(npts, 1);

  if (npts == 0)
    return;
  end

  pts = pts([end 1:end 1], :);

  for i=1:npts
    angle_prev = atan2(-(pts(i+1, 2) - pts(i, 2)), pts(i+1,1) - pts(i, 1));
    angle_next = atan2(-(pts(i+2, 2) - pts(i+1, 2)), pts(i+2,1) - pts(i+1, 1));
    angle = pi - (angle_next - angle_prev);

    if (angle > 2*pi)
      angle = angle - 2*pi;
    elseif (angle < 0)
      angle = angle + 2*pi;
    end

    interm = (pts(i+2, :) + pts(i, :)) / 2;
    greater_x = pts(1:end-1, 1) > interm(1, 1);

    indx = find(xor(greater_x(1:end-1), greater_x(2:end)));
    intersection = ((pts(indx+1, 2) - pts(indx, 2)) ./ (pts(indx+1, 1) - pts(indx, 1))) .* (interm(1, 1) - pts(indx, 1)) + pts(indx, 2);

    ninter = sum(intersection >= interm(1, 2), 1);

    if (mod(ninter, 2) == 0 & angle <= (pi - thresh) & angle >= thresh)
      conc(i) = true;
    end
  end

  if (sum(conc) == 1)
    indx = find(conc);
    half = round(npts/2);

    if (indx > half)
      indx = indx - half;
      if (indx < 1)
          indx = 1;
      end
      conc(indx) = true;
    else
      indx = indx + half;
      if (indx > npts)
          indx = npts;
      end
      conc(indx) = true;
    end
  elseif (sum(conc) == 0)
    conc([1 end]) = true;
  end

  return;
end
