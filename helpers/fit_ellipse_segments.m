function [ellipses, scores] = fit_ellipse_segments(pts, junctions, max_ratio, max_dist, max_score, max_overlap)

  nsegments = length(junctions);
  ellipses = NaN(nsegments, 5);
  npts = size(pts, 1);
  segments = cell(nsegments, 1);
  scores = Inf(nsegments, 1);

  for i=1:nsegments

    if (i == nsegments)
      index = [junctions(i):npts 1:junctions(1)];
    else
      index = [junctions(i):junctions(i+1)];
    end
    tmp = pts(index, :);

    if (size(tmp, 1) < 2 | (any(all(bsxfun(@eq, tmp(2:end, :), tmp(1,:)), 1), 2)) | all(diff(abs(diff(tmp)), [], 2) == 0))
      continue;
    end

    [ellipse, score] = fit_distance(tmp);
    
    segments{i} = tmp;
    ellipses(i,:) = ellipse;

    scores(i) = score;
  end

  [ellipses, segments, scores] = combine_ellipses(segments, ellipses, scores, max_ratio, max_dist, max_score, max_overlap);

  goods = (scores <= max_score & (ellipses(:, 4) ./ ellipses(:, 3)) >= max_ratio & ~any(isnan(ellipses),2));
  %ellipses = ellipses(scores <= max_score & (ellipses(:, 4) ./ ellipses(:, 3)) >= max_ratio, :);
  %ellipses = ellipses(~any(isnan(ellipses), 2), :);
  ellipses = ellipses(goods, :);
  scores = scores(goods);

  return;
end

function [ellipses, segments, scores] = combine_ellipses(segments, ellipses, scores, max_ratio, max_dist, max_score, max_overlap)

  nsegments = length(segments);
  improved = false;
  max_pow = max_dist^2;

  for i=1:nsegments
    if (isinf(scores(i)))
      continue
    end

    for j=i+1:nsegments
      if (isinf(scores(j)))
        continue
      end
      [ellipse, avg] = fit_distance([segments{i}; segments{j}]);

      % Ellipse selection
      if (avg > max_score | (ellipse(4) / ellipse(3)) < max_ratio | overlaps(ellipses, ellipse, [i,j]) > max_overlap)
        continue;

      % Case 1
      elseif (all(sum(bsxfun(@minus, ellipses([i j], 1:2), ellipse(1:2)).^2) > max_pow) | sum((ellipses(i, 1:2) - ellipses(j, 1:2)).^2) > 9*max_pow)
        continue;

      % Case 2
      % Case 3
      elseif ((all(ellipses([i j], 4) < max_dist) | all(abs(diff(ellipses([i j], 3:4))) < [1 0.5]*max_dist) | (abs(diff(ellipses([i j], 4)) < 0.05*max_dist))) | ...
             (avg < mean(scores([i j])) + std(scores([i j]))))
        segments{i} = [segments{i}; segments{j}];
        segments{j} = [];
        scores(i) = avg;
        scores(j) = Inf;
        ellipses(j, :) = NaN;
        ellipses(i, :) = ellipse;
        improved = true;
      end
    end
  end

  if (improved)
    [ellipses, segments, scores] = combine_ellipses(segments, ellipses, scores, max_ratio, max_dist, max_score, max_overlap);
  end

  return;
end

function [overlap] = overlaps(ellipses, ellipse, indxs)

  if (numel(indxs) == 1 & indxs == -1)
    overlap = zeros(size(ellipses, 1), 1);
  else
    overlap = 0;
  end
  curr_area = prod(ellipse(3:4))*pi;

  for i=1:size(ellipses, 1)
    if (all(i ~= indxs) & ~isnan(ellipses(i, 1)))
      tmp = ellipse_ellipse_area_mex(ellipse(5), ellipse(3), ellipse(4), ellipse(1), ellipse(2), ellipses(i,5), ellipses(i,3), ellipses(i,4), ellipses(i,1), ellipses(i,2));

      if (numel(overlap) == 1)
        tmp = tmp / curr_area;

        if (tmp > overlap)
          overlap = tmp;
        end
      else
        overlap(i) = tmp;
      end
    end
  end

  return;
end

function [ellipse, avg] = fit_distance(pts)

  ellipse = NaN(1, 5);
  avg = Inf;

  if (isempty(pts))
    return;
  end

  [c, a, o] = fit_ellipse(pts);
  ell_pts = carth2elliptic(pts, c, a, o, 'radial');

  dist = abs(ell_pts(:,2) - 1);

  thresh = 4*std(dist);
  pts = pts(dist < thresh, :);

  if (isempty(pts))
    dist = [];

    return;
  end

  [c, a, o] = fit_ellipse(pts);
  if (any(isnan(c)))
    return;
  end
  ellipse = [c.' a.' o];
  ell_pts = carth2elliptic(pts, c, a, o, 'radial');

  dist = abs(ell_pts(:,2) - 1);
  avg = mean(dist);

  return;
end
