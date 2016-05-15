function [mymovie, all_estim] = detect_embryos(mymovie, estim_only, opts)

  min_egg_size = [7 5] um;


  if (nargin == 2)
    opts = estim_only;
    estim_only = false;
  end

  type = opts.segmentation_type;

  if (isstruct(mymovie))
    switch type
      case 'dic'
        [nframes, imgsize] = size_data(mymovie.dic);
      case 'all'
        if (isfield(mymovie, 'eggshell') & ~isempty(mymovie.eggshell))
          opts.segmentation_type = 'markers';
          mymovie = split_cells(mymovie, estim_only, opts);
        end
        opts.segmentation_type = 'dic';
        mymovie = split_cells(mymovie, estim_only, opts);

        return;
      case 'markers'
        if (~isfield(mymovie, 'eggshell') | isempty(mymovie.eggshell))
          opts.segmentation_type = 'dic';
          mymovie = split_cells(mymovie, estim_only, opts);

          return;
        end
        type = 'eggshell';
        [nframes, imgsize] = size_data(mymovie.eggshell);
      case 'data'
        [nframes, imgsize] = size_data(mymovie.data);
      otherwise
        error 'None of the expected field are present in ''mymovie''';
    end
  else
    [h,w,nframes] = size(imgs);
    imgsize = [h, w];
  end

  rad_thresh = max(egg_min_size);
  rad_thresh = rad_thresh ./ [0.5 4 8 10 12];
  area_thresh = prod(egg_min_size + rad_thresh(3))*pi;

  opts.split_parameters.max_distance = opts.split_parameters.max_distance / opts.pixel_size;

  npixels = max(imgsize);
  size10 = round(npixels/10);
  size75 = round(npixels/75);
  size100 = round(prod(imgsize) / 100);
  size150 = round(npixels/150);
  size250 = round(npixels/250);
  size200 = round(npixels/200);

  all_ellipses = cell(nframes, 1);
  all_estim = cell(nframes, 1);

  max_nellipses = 0;
  ndata = 0;

  for n = 1:nframes
    %nimg = randi(nframes, 1);
    nimg = n;
    %nimg = 8

    if (isstruct(mymovie))
      
      switch type
        case 'dic'
          if (opts.verbosity == 3)
            img = double(load_data(mymovie.dic, nimg));
            figure;imshow(imnorm(img));
            img = imadm_mex(img);
            figure;imshow(imnorm(img));
            thresh = graythresh(img);
            img = (img > thresh*0.5*(max(img(:))) );
            figure;imshow(img);
          end
          img = double(load_data(mymovie.dic, nimg));
          img = imadm_mex(img);
          thresh = graythresh(img);
          img = (img > thresh*0.5*(max(img(:))) );
        case 'eggshell'
          img = double(load_data(mymovie.eggshell, nimg));

          noise_params = estimate_noise(img);

          img = gaussian_mex(img, size250);
          img = median_mex(img, size200, 3);

          img = (img < noise_params(1) + noise_params(2));

          %thresh = graythresh(img);          
          %if (thresh == 0)
          %  warning('Might need to normalize the image to detect the embryos');
          %end
          %img = (img > thresh*(max(img(:))));
        case 'data'
          %img = imnorm(double(load_data(mymovie.data, nimg)));
          img = (double(load_data(mymovie.data, nimg)));
          noise_params = estimate_noise(img);
          vals = range(img(:));

          noise_thresh = min(round(vals / (15*noise_params(2))) + 1, 10);

          img = gaussian_mex(img, size250);
          img = median_mex(img, size200, 3);

          %img = (img > noise_params(1) + 10*noise_params(2));
          img = (img > noise_params(1) + noise_thresh*noise_params(2));
          %thresh = graythresh(img);
          %img = (img > thresh*(max(img(:))) );
        otherwise
          error 'Not implemented yet'
      end
    else
      img = mymovie(:, :, nimg);
    end

    %imagesc(img);title(num2str(n));drawnow;
    %keyboard

    img = padarray(img, [size10 size10]);

    img = imdilate(img, strel('disk', size150));
    img = imfill(img, 'holes');
    img = bwareaopen(img, size100);
    img = imdilate(img, strel('disk', size75));
    img = imfill(img, 'holes');
    img = imerode(img, strel('disk', size75 + size150));

    img =  img((size10+1):(size10+imgsize(1)),(size10+1):(size10+imgsize(2)));

    if (opts.verbosity == 3)
      figure;imshow(img);
    end

    if(~any(img) | (sum(img(:)) / prod(imgsize) > 0.9))
      continue;
    end

    estim = bwboundaries(img, 8, 'noholes');
    if (isempty(estim))
      continue;
    end

    if (estim_only)
      all_estim{nimg} = estim{1};
      continue;
    end

    for i=1:length(estim)
      tmp_estim = estim{i};
      tmp_estim = tmp_estim(:,[2 1]);

      [pac, indxs] = impac(tmp_estim);

      %borders = (any(tmp_estim == 2 | bsxfun(@eq, tmp_estim, imgsize([2 1])-1), 2));
      borders = (any(tmp_estim == 2 | bsxfun(@eq, tmp_estim, imgsize([2 1])-1), 2));
      border_indx = find(xor(borders, borders([2:end 1])));

      concaves = compute_concavity(pac, opts.split_parameters.angle_thresh);

      [indxs, indx_indx] = sort([indxs; border_indx]);
      concaves = [concaves; true(size(border_indx))];
      concaves = concaves(indx_indx);

      ellipses = fit_segments(tmp_estim, indxs(concaves), borders, opts.split_parameters.max_ratio, opts.split_parameters.max_distance, opts.split_parameters.max_score, opts.split_parameters.max_overlap);

      if (opts.verbosity == 3)
        figure
        imshow(img);
        hold on
        scatter(tmp_estim(borders, 1), tmp_estim(borders, 2), 'g');
        scatter(tmp_estim(indxs, 1), tmp_estim(indxs, 2), 'r');
        scatter(tmp_estim(indxs(concaves), 1), tmp_estim(indxs(concaves), 2), 'y');

        hold on;
        for j = 1:size(ellipses, 1)
          draw_ellipse(ellipses(j, 1:2), ellipses(j, 3:4), ellipses(j, 5));
        end
      end

      if (i == 1)
        all_ellipses{nimg} = ellipses;
        all_estim{nimg} = tmp_estim;
      else
        all_ellipses{nimg} = [all_ellipses{nimg}; ellipses];
        all_estim{nimg} = [all_estim{nimg}; [NaN, NaN]; tmp_estim];
      end
    end

    if (ndata == 0)
      ndata = size(all_ellipses{nimg}, 2);
    end
    nellipses = size(all_ellipses{nimg}, 1);
    if (nellipses > max_nellipses)
      max_nellipses = nellipses;
    end

    if (i > 1)
      all_estim{nimg} = [all_estim{nimg}; [NaN, NaN]];
    end
  end

  if (nframes == 1)
    all_ellipses = all_ellipses{1};
    all_estim = all_estim{1};

  elseif (max_nellipses > 1)
    display('There seems to be several cells in the images, double checking if everything is fine !');

    dist_thresh = (opts.split_parameters.max_distance / 3)^2;
    start_index = nframes;
    for i=1:nframes
      if (~isempty(all_ellipses{i}))
        real_ell = NaN(size(all_ellipses{i},1), ndata, nframes);
        real_ell(:,:,1) = all_ellipses{i};
        start_index = i;

        break;
      end
    end

    for i=start_index+1:nframes
      tmp_ell = all_ellipses{i};
      for k=1:size(tmp_ell, 1)
        found = false;
        mpos = mymean(real_ell(:,1:2,:), 3);
        for j=1:size(mpos, 1)
          if (sum((mpos(j, :) - tmp_ell(k, 1:2)).^2) < dist_thresh)
            found = true;
            real_ell(j, :, i) = tmp_ell(k, :);
            break;
          end
        end

        if (~found)
          real_ell(end+1,:,:) = NaN;
          real_ell(end, :, 1) = tmp_ell(k,:);
        end
      end
    end

    goods = sum(~isnan(real_ell),3);
    goods = (goods(:,1) > nframes/3);

    avg_ell = NaN(0,ndata);

    for i=1:size(real_ell, 1)
      if (goods(i))
        tmp = real_ell(i, :, :);
        tmp = reshape(tmp(~isnan(tmp)), 5, []);
        avg_ell(end+1,:) = median(tmp, 2);

        %draw_ellipse(avg_ell(end, 1:2), avg_ell(end, 3:4), avg_ell(end, 5), 'r');
      end
    end
    avg_area = prod(avg_ell(:,3:4), 2)*pi;
    mean_area = mean(avg_area);
    goods = (avg_area <= avg_area * opts.split_parameters.max_area_diff & ...
             avg_area >= avg_area / opts.split_parameters.max_area_diff);

    avg_ell = avg_ell(goods, :);
    ntargets = size(avg_ell, 1);

    for nimg = 1:nframes
      
      estim = all_estim{nimg};
      
      if (isempty(estim))
        continue;
      end
          
      [pac, indxs] = impac(estim);
      concaves = compute_concavity(pac, opts.split_parameters.angle_thresh);

      sides = (any(estim == 2 | bsxfun(@eq, estim, imgsize-1), 2));
      borders = (xor(sides, sides([2:end 1])));
      sides(borders) = false;

      borders = borders | isnan(estim(:,1));
      borders(indxs(concaves)) = true;
      
      labels = cumsum(double(borders), 1);

      estim = estim(~sides, :);
      labels = labels(~sides);

      nlabel = hist(labels, labels(end));
      label_indx = [1:length(nlabel)];

      too_few = ismember(labels, label_indx(nlabel < 10));
      estim = estim(~too_few, :);
      labels = labels(~too_few);

      indxs = unique(labels);

      ellipses = NaN(ntargets, 5);
      scores = Inf(ntargets, 1);

      for i = 1:ntargets
        %current = ismember(labels, indxs);
        %tmp_pts = estim(current, :);
        ell_pts = carth2elliptic(estim, avg_ell(i, 1:2), avg_ell(i, 3:4), avg_ell(i, 5), 'radial');
        dist = abs(ell_pts(:,2) - 1);

        valids = (dist <= 2*opts.split_parameters.max_score);
        npts = mymean(double(valids), 1, labels);

        fit_indxs = indxs(npts > 0.333);
        current = ismember(labels, fit_indxs);
        [ellipses(i, :)] = fit_distance(estim(current, :));
        scores(i) = overlaps(avg_ell(i,:), ellipses(i,:), NaN); 
      end

      goods = (scores > 0.75 & (ellipses(:, 4) ./ ellipses(:, 3)) >= opts.split_parameters.max_ratio);
      ellipses(~goods, :) = avg_ell(~goods, :);

       % for j = 1:size(ellipses, 1)
       %   draw_ellipse(ellipses(j, 1:2), ellipses(j, 3:4), ellipses(j, 5));
       % end

        all_ellipses{nimg} = ellipses;
    end
  end

  if (isstruct(mymovie))
    neighbors = get_struct('reference');

    for i=1:nframes
      ellipse = all_ellipses{i};
      
      if (isempty(ellipse))
        neighbors.centers = NaN(2,1);
        neighbors.axes_length = NaN(2,1);
        neighbors.orientations = NaN;
      else
        [neighbors.centers, neighbors.axes_length, neighbors.orientations] = deal(ellipse(:, 1:2).', ellipse(:, 3:4).', ellipse(:, 5).');

        dist = [bsxfun(@minus, ellipse(:, [1 2]), imgsize([2 1])).^2 ellipse(:, [1 2]).^2];
        dist = min(dist, [], 2);
        [junk, indx] = max(dist);
        neighbors.index = indx;
      end
      
      if (isempty(mymovie.(type).neighbors))
        mymovie.(type).neighbors = neighbors;
      else
        mymovie.(type).neighbors(i) = neighbors;
      end
    end
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
