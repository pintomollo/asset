function img = mask_neighbors(img, neighbors, opts)
  
  %keyboard
  
  nneigh = size(neighbors.centers, 2);
  if (nneigh == 0|all(isnan(neighbors.centers(:))))
    return;
  end

  alpha = 0.075;
  beta =  2.75 / opts.pixel_size;
  %gamma = params.gamma;

  growth_size = ceil(0.75 / opts.pixel_size);
  filt = strel('disk', growth_size);

  mask = false(size(img));
  for i=1:nneigh
    path = draw_ellipse(neighbors.centers(:, i), neighbors.axes_length(:, i), neighbors.orientations(i));
    mask = mask | roipoly(mask, path(:,1), path(:,2));
  end

  repl_val = median(img(:));

  mask = imdilate(mask, filt);
  dist = bwdist(~mask);
  dist = 1 ./ (1 + exp(-alpha * (dist - beta)));
  %weight(weight > 1 - gamma) = Inf;
  dist(~mask) = 0;

  dist = imnorm(dist);
  img = img .* (1-dist) + repl_val .* dist;

  return;
end
