function spots = tracking(spots, track_opts, opts)

  % Input checking
  if (nargin < 2)
    opts = get_struct('ASSET');
    track_opts = opts.spot_tracking;
  elseif (nargin == 2)
    if (isfield(track_opts, 'linking_function'))
      opts = get_struct('ASSET');
    else
      opts = track_opts;
      track_opts = opts.spot_tracking;
    end
  end

  nframes = length(spots);

  % Initialize the output variable
  links = cell(nframes, 1);

  spot_max_movement = track_opts.frame_displacement;
  frame_weight = track_opts.linking_function;
  max_frames = track_opts.frame_window;

  % If the size of the pixels has been set, we can compute the actual spot size
  if (opts.pixel_size > 0)

    % The maximal size of the spot in pixels
    spot_max_movement = (spot_max_movement / opts.pixel_size);
  else

    % Try to compute the pixel size
    opts = set_pixel_size(opts);

    % If it worked, we can now compute the size in pixels
    if (opts.pixel_size > 0)
      spot_max_movement = (spot_max_movement / opts.pixel_size);
    end
  end

  pts = [];
  npts = 0;
  ndim = 0;

  all_pts = cell(nframes, 1);
  all_assign = [];
  prev_max = NaN;

  for i=1:nframes

    prev_pts = pts;
    prev_npts = npts;

    pts = [spots(i).carth/opts.pixel_size spots(i).properties];
    pts = pts(~any(isnan(pts), 2), :);
    all_pts{i} = pts;
    [npts, tmp] = size(pts);

    if (ndim == 0)
      ndim = tmp;
    end

    if (prev_npts > 0 & npts > 0)

      dist = Inf(npts + prev_npts);

      mutual_dist = frame_weight(prev_pts, pts);
      mutual_dist(mutual_dist > spot_max_movement) = Inf;

      %%%%%%%%% Not updated properly as the article states...
      prev_max = spot_max_movement;

      ends = eye(max(npts,prev_npts))*prev_max;
      ends(ends==0) = Inf;
      trans_dist = (mutual_dist.' < Inf) * min(mutual_dist(:));
      trans_dist(trans_dist == 0) = Inf;

      dist(1:prev_npts,1:npts) = mutual_dist;
      dist(1:prev_npts,npts+1:end) = ends(1:prev_npts,1:prev_npts);
      dist(prev_npts+1:end,1:npts) = ends(1:npts,1:npts);
      dist(prev_npts+1:end,npts+1:end) = trans_dist;

      [assign, cost] = munkres(dist);
      assign_dist = dist(sub2ind(size(dist),[1:prev_npts+npts],assign));
      assign_dist = assign_dist(:);
      good_indx = (assign(1:prev_npts) <= npts);
      all_assign = [all_assign; (assign_dist(good_indx))];

      [junk tmp] = sort(assign);

      valids = (tmp <= prev_npts);
      valids(npts+1:end) = false;
      links{i} = [junk(valids).' tmp(valids).' (i-ones(sum(valids), 1))];
    end

    if (isempty(links{i}))
      links{i} = NaN(0, 3);
    end
  end

  if (max_frames == 0)
    return;
  end

  avg_movement = mean(all_assign);

  starts = zeros(0, ndim+2);
  ends = zeros(0, ndim+2);
  interm = zeros(0, ndim+2);

  joining_weight = track_opts.joining_function;
  splitting_weight = track_opts.splitting_function;
  prev_starts = [];

  for i=nframes:-1:2

    indx_interm = links{i}(:,1);
    indx_starts = [];

    nstarts = size(all_pts{i},1);
    indx_starts = setdiff([1:nstarts], indx_interm);
    nstarts = length(indx_starts);

    if (nstarts>0)
      starts(end+1:end+nstarts,:) = [all_pts{i}(indx_starts,:) indx_starts(:) ones(nstarts,1)*i];
    end

    nends = size(all_pts{i-1},1);
    indx_ends = setdiff([1:nends], links{i}(:,2));
    nends = length(indx_ends);
    if (nends>0)
      ends(end+1:end+nends,:) = [all_pts{i-1}(indx_ends,:) indx_ends(:) ones(nends,1)*i-1];
    end

    if (~isempty(all_pts{i}))

      join_dist = frame_weight(all_pts{i-1}(indx_ends,:), all_pts{i}(indx_interm,:));

      if (~isempty(prev_starts))
        split_dist = frame_weight(all_pts{i+1}(prev_starts, :), all_pts{i}(indx_interm, :));
        join_dist = [join_dist; split_dist];
      end

      indx_interm = indx_interm(any(join_dist < spot_max_movement, 1));
      ninterm = length(indx_interm);

      if (ninterm>0)
        interm(end+1:end+ninterm,:) = [all_pts{i}(indx_interm,:) indx_interm(:) ones(ninterm,1)*i];
      end
    end

    prev_starts = indx_starts;
  end

  nstarts = size(starts, 1);
  nends = size(ends, 1);
  ninterm = size(interm, 1);

  dist = Inf(nstarts + nends + 2*ninterm);

  closing_func = track_opts.gap_function;

  mutual_dist = closing_func(ends, starts);

  frame_indx = -bsxfun(@minus,ends(:,end),starts(:,end).');

  mutual_dist(frame_indx < 1 | frame_indx > max_frames | mutual_dist > spot_max_movement) = Inf;

  alt_cost = prctile(mutual_dist(~isinf(mutual_dist)), 90);
  alt_dist = eye(max(nends, nstarts))*alt_cost;
  alt_dist(alt_dist==0) = Inf;

  [merge_dist, merge_weight, alt_weight] = joining_weight(ends, interm, all_pts, links);
  frame_indx = -bsxfun(@minus, ends(:,end), interm(:,end).');
  merge_weight = merge_dist .* merge_weight;
  merge_weight(merge_dist > spot_max_movement | frame_indx ~= 1) = Inf;
  alt_weight = avg_movement * alt_weight;

  [split_dist, split_weight, alt_split_weight] = splitting_weight(interm, starts, all_pts, links);
  frame_indx = -bsxfun(@minus, interm(:,end), starts(:,end).');
  split_weight = split_dist .* split_weight;
  split_weight(split_dist > spot_max_movement | frame_indx ~= 1) = Inf;
  alt_split_weight = avg_movement * alt_split_weight;

  % Note that end-end merging and start-start splitting is not allowed by this
  % algorithm, which might make sense...

  dist(1:nends,1:nstarts) = mutual_dist;
  dist(1:nends,nstarts+1:nstarts+ninterm) = merge_weight;
  dist(nends+1:nends+ninterm,1:nstarts) = split_weight;

  dist(1:nends,nstarts+ninterm+1:end-ninterm) = alt_dist(1:nends,1:nends);
  dist(nends+ninterm+1:end-ninterm,1:nstarts) = alt_dist(1:nstarts,1:nstarts);

  dist(nends+1:nends+ninterm,end-ninterm+1:end) = alt_split_weight;
  dist(end-ninterm+1:end,nstarts+1:nstarts+ninterm) = alt_weight;

  trans_dist = (dist(1:nends+ninterm, 1:nstarts+ninterm).' < Inf) * min(mutual_dist(:));
  trans_dist(trans_dist == 0) = Inf;

  dist(nends+ninterm+1:end,nstarts+ninterm+1:end) = trans_dist;

  [assign, cost] = munkres(dist);

  for i=1:nends+ninterm
    if (assign(i) <= nstarts)
      target = starts(assign(i), :);
    elseif (assign(i) < nstarts + ninterm)
      target = interm(assign(i) - nstarts, :);
    else
      continue;
    end

    if (i <= nends)
      reference = [target(end-1) ends(i,end-1:end)];
    else
      reference = [target(end-1) interm(i-nends,end-1:end)];
    end
    links{target(end)} = [links{target(end)}; reference];
  end

  for i=1:nframes
    spots(i).cluster = links{i};
  end

  return
end
