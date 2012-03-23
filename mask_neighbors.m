function img = mask_neighbors(img, center, axes_length, orientation, neighbors, opts)
  
  nneigh = size(neighbors.centers, 2);
  if (nneigh <= 0 | all(isnan(neighbors.centers(:))))
    return;
  end

  alpha = 0.075;
  beta =  1 / opts.pixel_size;
  %gamma = params.gamma;

  growth_size = ceil(1.25 / opts.pixel_size);
  filt = strel('disk', growth_size);

  if (opts.verbosity == 3)
    figure;imshow(img);
    hold on;
  end

  mask = false(size(img));
  for i=1:nneigh
    if (i ~= neighbors.index)
      path = draw_ellipse(neighbors.centers(:, i), neighbors.axes_length(:, i) + growth_size, neighbors.orientations(i));
      ell_path = carth2elliptic(path, center, axes_length, orientation, 'radial');
      path = path(ell_path(:,2) >= 1, :);
      mask = mask | roipoly(mask, path(:,1), path(:,2));

      if (opts.verbosity == 3)
        plot(path(:,1), path(:,2));
      end
    end
  end

  if (opts.verbosity == 3)
    figure;imshow(mask)
  end

  repl_val = median(img(:));

  %mask = imdilate(mask, filt);
  dist = double(bwdist(~mask));
  dist = 1 ./ (1 + exp(-alpha * (dist - beta)));
  %weight(weight > 1 - gamma) = Inf;
  dist(~mask) = 0;

  dist = imnorm(dist);
  
  if (opts.verbosity == 3)
    figure;imshow(dist)
  end

  img = img .* (1-dist) + repl_val .* dist;

  if (opts.verbosity == 3)
    figure;imshow(img)
  end

  return;
end
