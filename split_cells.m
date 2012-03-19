function [all_ellipses, all_estim] = split_cells(imgs, estim_only, opts, params)

  if (nargin == 2)
    opts = estim_only;
    estim_only = false;
  end

  [h,w,nframes] = size(imgs);
  imgsize = [w, h];

  if (nargin == 4)
    [max_ratio, angle_thresh, max_score, max_overlap, max_dist] = deal(params(1), params(2), params(3), params(4), params(5));
  else
    max_ratio = 0.5;
    angle_thresh = 0.2;
    max_score = 0.05;
    max_overlap = 0.4;

    max_dist = 15;
  end

  max_dist = max_dist / opts.pixel_size;

  npixels = max(imgsize);
  size10 = round(npixels/10);
  size75 = round(npixels/75);
  size100 = round(h*w / 100);
  size150 = round(npixels/150);

  all_ellipses = cell(nframes, 1);
  all_estim = cell(nframes, 1);

  for n = 1:nframes
    %nimg = randi(nframes, 1);
    nimg = n;
    %nimg = 23

    img = imgs(:, :, nimg);

    %imagesc(img);hold on;
    img = padarray(img, [size10 size10]);

    img = imdilate(img, strel('disk', size150));
    img = imfill(img, 'holes');
    img = bwareaopen(img, size100);
    img = imdilate(img, strel('disk', size75));
    img = imfill(img, 'holes');
    img = imerode(img, strel('disk', size75 + size150));

    img =  img((size10+1):(size10+h),(size10+1):(size10+w));

    if(~any(img) | (sum(img(:)) / prod(imgsize) > 0.9))
      continue;
    end

    estim = bwboundaries(img, 8, 'noholes');
    if (isempty(estim))
      continue;
    end

    if (estim_only)
      all_estim{n} = estim{1};
      continue;
    end

    for i=1:length(estim)
      tmp_estim = estim{i};
      tmp_estim = tmp_estim(:,[2 1]);

      [pac, indxs] = impac(tmp_estim);

      borders = (any(tmp_estim == 2 | bsxfun(@eq, tmp_estim, imgsize-1), 2));
      border_indx = find(xor(borders, borders([2:end 1])));

      concaves = compute_concavity(pac, angle_thresh);

      [indxs, indx_indx] = sort([indxs; border_indx]);
      concaves = [concaves; true(size(border_indx))];
      concaves = concaves(indx_indx);

      ellipses = fit_segments(tmp_estim, indxs(concaves), borders, max_ratio, max_dist, max_score, max_overlap);

      if (opts.verbosity == 3)
        figure
        imshow(img);
        hold on
        scatter(tmp_estim(borders, 1), tmp_estim(borders, 2), 'g');
        scatter(tmp_estim(indxs, 1), tmp_estim(indxs, 2), 'r');
        scatter(tmp_estim(indxs(concaves), 1), tmp_estim(indxs(concaves), 2), 'y');

        %hold on;
        for j = 1:size(ellipses, 1)
          draw_ellipse(ellipses(j, 1:2), ellipses(j, 3:4), ellipses(j, 5));
        end
      end

      if (i == 1)
        all_ellipses{n} = ellipses;
        all_estim{n} = tmp_estim;
      else
        all_ellipses{n} = [all_ellipses{n}; ellipses];
        all_estim{n} = [all_estim{n}; [NaN, NaN]; tmp_estim];
      end
    end
    if (i > 1)
      all_estim{n} = [all_estim{n}; [NaN, NaN]];
    end
  end

  if (nframes == 1)
    all_ellipses = all_ellipses{1};
    all_estim = all_estim{1};
  end

  if (nargout == 1)
    all_estim = [];
  end

  return;
end

function ellipses = fit_segments(pts, junctions, is_border, max_ratio, max_dist, max_score, max_overlap)

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
    tmp_border = is_border(index);
    tmp = tmp(~tmp_border, :);

    if (size(tmp, 1) < 2 | (any(all(bsxfun(@eq, tmp(2:end, :), tmp(1,:)), 1), 2)) | all(diff(abs(diff(tmp)), [], 2) == 0))
      continue;
    end

    [ellipse, score] = fit_distance(tmp);
    
    segments{i} = tmp;
    ellipses(i,:) = ellipse;
    scores(i) = score;
  end

  [ellipses, segments, scores] = combine_ellipses(segments, ellipses, scores, max_ratio, max_dist, max_score, max_overlap);

  ellipses = ellipses(scores <= max_score & (ellipses(:, 4) ./ ellipses(:, 3)) >= max_ratio, :);
  ellipses = ellipses(~any(isnan(ellipses), 2), :);

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

  overlap = 0;
  curr_area = prod(ellipse(3:4))*pi;

  for i=1:size(ellipses, 1)
    if (all(i ~= indxs) & ~isnan(ellipses(i, 1)))
      tmp = ellipse_ellipse_area_mex(ellipse(5), ellipse(3), ellipse(4), ellipse(1), ellipse(2), ellipses(i,5), ellipses(i,3), ellipses(i,4), ellipses(i,1), ellipses(i,2));
      tmp = tmp / curr_area;

      if (tmp > overlap)
        overlap = tmp;
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
      conc(indx - half) = true;
    else
      conc(indx + half) = true;
    end
  elseif (sum(conc) == 0)
    conc([1 end]) = true;
  end

  return;
end
