function [path, emissions] = remove_polar_body(img, path, params, scoring_func, weights, opts, emissions)

  [height, width] = size(img);
  repeats = ceil(height / 4);

  vert_indexes = [1:height]';

  opts.force_circularity = false;
  opts.dp_method = 'normal';
  weights_path = weights.path;
  weights.alpha = 1;

  max_dist = height / 10;
  security = round(width / 100);
  coef = 4;

  dx = diff(path);

  if (opts.verbosity == 3)
    old_path = path;

    figure;plot(dx)
    hold on;
    plot([1 1; length(dx) length(dx)], [1 -1; 1 -1]*(mean(dx) + coef*std(dx)))
  end

  indx = find(abs(dx) > mean(dx) + coef*std(dx));

  if (numel(indx) < 2)
    return;
  end

  indx = fuse_indexes(indx, dx);

  %indx = [indx; indx(1)+height]; 

  %dist = diff(indx);
  %cands = find(dist < max_dist);

  for i=1:length(indx)
    for j=1:length(indx)
      if (i ~= j & dx(indx(i)) > 0 & dx(indx(j)) < 0)
        if (i > j)
          around = true;
        else
          around = false;
        end

        if (indx(j) - indx(i) + height*around < max_dist)

          myrange = [indx(i) indx(j)];
          normals = find(abs(dx) < mean(dx) + std(dx));
          start = (myrange(1) - normals);
          start = start(start > 0);
          start = -min(start);

          ends = (normals - myrange(2));
          ends = ends(ends > 0);
          ends = min(ends);

          myrange = myrange + [start ends];
          if (around)
            full_range = [myrange(1):height 1:myrange(2)];
          else
            full_range = [myrange(1):myrange(2)];
          end

          rad_range = [min(path(myrange(1)), path(myrange(2))) max(path(myrange(1)), path(myrange(2)))] + [-security security];
          rad_range(rad_range < 0) = 1;
          rad_range(rad_range > width) = width;
          
          subimg = img(full_range, rad_range(1):rad_range(2));
          params.init = round(path(myrange(1)) - rad_range(1) + 1);
          params.final = round(path(myrange(2)) - rad_range(1) + 1);

          if (~isempty(weights_path))
            weights.path = weights_path(full_range,:) - rad_range(1) + 1;
          end

          if (nargin < 7)
            cortex_path = dynamic_programming(subimg, params, scoring_func, weights, opts);

            if (opts.verbosity == 3)
              figure;
              imshow(subimg);
              hold on;plot(cortex_path, [1:length(cortex_path)])
            end
          else
            [cortex_path, new_emission, junk] = dynamic_programming(subimg, params, scoring_func, weights, opts);
            emissions(full_range, rad_range(1):rad_range(2)) = new_emission;
          end
          
          path(full_range) = cortex_path + rad_range(1) - 1;      
        end
      end
    end
  end

  if (opts.verbosity == 3)
    figure;imshow(img);
    hold on;plot(old_path, vert_indexes);
    plot(path, vert_indexes, 'r');
  end

  return;
end

function indx = fuse_indexes(indx, values)

  prev_val = values(indx(end));
  prev_indx = indx(end);
  for i=length(indx)-1:-1:1
    if (prev_indx - indx(i) < 5 & sign(prev_val) == sign(values(indx(i))))
      prev_indx = indx(i);

      if (sign(prev_val) > 0)
        indx(i+1) = [];
      else
        indx(i) = [];
      end
    else
      prev_indx = indx(i);
      prev_val = values(indx(i));
    end
  end

  return;
end
