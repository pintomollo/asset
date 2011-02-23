function img = erase_egg(img, outer, inner)

  [nrows,npts] = size(img);

  if (length(outer) ~= nrows | length(inner) ~= nrows)
    return;
  end

  outer = round(2*outer - inner);

  npixels = npts - outer;
  starts = round(inner);
  ends = starts + npixels;

  if (any(isnan(starts) | isinf(starts)))
    return;
  end

  for i=1:nrows
    img(i,starts(i,1):ends(i,1)) = img(i, outer(i,1):end);
  end

  return;
end
