function img = binning(img, factor)

  if (factor <= 1)
    return;
  end

  img = img(1:2:end, :) + img (2:2:end, :);
  img = img(:, 1:2:end) + img (:, 2:2:end);

  factor = factor / 2;

  img = binning(img, factor);

  return;
end
