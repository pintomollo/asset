function [links, spots] = track_spots(spots, opts)

  % Input checking
  if (nargin < 2)
    opts = get_struct('ASSET');
  end

  if (isstruct(spots))
    mymovie = spots;
    if (isfield(mymovie, 'experiment'))
      tmp_struct = mymovie.data.spots;
    else
      tmp_struct = mymovie;
    end

    nframes = length(tmp_struct);
    spots = cell(nframes, 1);

    for i=1:nframes
      spots{i} = [tmp_struct(i).carth tmp_struct(i).properties];
    end
  else
    mymovie = [];
    nframes = length(spots);
  end

  % Initialize the output variable
  links = cell(nframes, 1);

  if (isfield(opts, 'pixel_size'))
    opts = set_pixel_size(opts);

    pixel_size = opts.pixel_size;
  elseif (isfield(opts, 'spot_tracking') & isfield(opts.spot_tracking, 'pixel_size'))
    opts.spot_tracking = set_pixel_size(opts.spot_tracking);

    pixel_size = opts.spot_tracking.pixel_size;
  else
    pixel_size = 1;
  end

  if (pixel_size <= 0)
    pixel_size = 1;
  end

  if (isfield(opts, 'spot_tracking'))
    opts = opts.spot_tracking;
  end

  spot_max_movement = opts.frame_displacement;
  frame_weight = opts.linking_function;
  max_frames = opts.frame_window;

  % The maximal size of the spot in pixels
  spot_max_movement = (spot_max_movement / pixel_size);

  pts = [];
  npts = 0;
  ndim = 0;

  all_assign = [];
  prev_max = NaN;

  do_display = (isfield(opts, 'verbosity') && opts.verbosity > 1);
  if (do_display)
    cpb = ConsoleProgressBar();
    cpb.setLeftMargin(4);   % progress bar left margin
    cpb.setTopMargin(1);    % rows margin
    cpb.setLength(100);      % progress bar length: [.....]
    cpb.setMinimum(0);
    cpb.setMaximum(nframes);

    cpb.setElapsedTimeVisible(1);
    cpb.setRemainedTimeVisible(1);
    cpb.setElapsedTimePosition('left');
    cpb.setRemainedTimePosition('right');

    cpb.start();
  end

  nprops = 0;
  for i=1:nframes

    prev_pts = pts;
    prev_npts = npts;

    pts = spots{i};
    [npts, tmp] = size(pts);

    if (npts > 0 && nprops == 0)
      nprops = tmp - 2;
    end

    if (ndim == 0)
      ndim = tmp;
    end
    
    if (prev_npts > 0 & npts > 0)

      dist = Inf(npts + prev_npts);
      
      [mutual_dist, weight] = frame_weight(prev_pts, pts);
      mutual_dist(mutual_dist > spot_max_movement) = Inf;
      mutual_dist = mutual_dist .* weight;
      prev_max = spot_max_movement;

      ends = eye(max(npts,prev_npts))*prev_max;
      ends(ends==0) = Inf;
      min_dist = min(mutual_dist(:));
      trans_dist = Inf(npts,prev_npts);
      trans_dist(mutual_dist.' < Inf) = min_dist;

      dist(1:prev_npts,1:npts) = mutual_dist;
      dist(1:prev_npts,npts+1:end) = ends(1:prev_npts,1:prev_npts);
      dist(prev_npts+1:end,1:npts) = ends(1:npts,1:npts);
      dist(prev_npts+1:end,npts+1:end) = trans_dist;

      %[assign, cost] = munkres(dist);
      %[assign, cost] = assignmentoptimal(dist);
      [assign, cost] = lapjv_fast(dist);

      assign = assign(1:prev_npts);
      indxs = 1:length(assign);
      good_indx = (assign <= npts);

      assign = assign(good_indx);
      indxs = indxs(good_indx);

      assign_dist = dist(sub2ind(size(dist),indxs,assign));
      %assign_dist = assign_dist(:);
      %good_indx = (assign(1:prev_npts) <= npts);
      all_assign = [all_assign; assign_dist(:)];
      
      [assign, perms] = sort(assign(:));
      %[junk tmp] = sort(assign);

      %valids = (tmp <= prev_npts);
      %valids(npts+1:end) = false;
      %links{i} = [junk(valids) tmp(valids) (i-ones(sum(valids), 1))];
      links{i} = [assign indxs(perms).' (i-ones(length(assign), 1))];
    end

    if (isempty(links{i}))
      links{i} = NaN(0, 3);
    end

    if (do_display)
      text = sprintf('Progress: %d/%d', i, nframes);
      cpb.setValue(i);
      cpb.setText(text);
    end
  end

  if (do_display)
    cpb.stop();
  end

  if (max_frames == 0)
    return;
  end

  avg_movement = mean(all_assign);

  starts = zeros(0, ndim+2);
  ends = zeros(0, ndim+2);
  interm = zeros(0, ndim+2);

  closing_weight = opts.gap_function;
  joining_weight = opts.joining_function;
  splitting_weight = opts.splitting_function;

  tracking_options = ~[isempty(closing_weight) isempty(joining_weight) isempty(splitting_weight)];
  prev_starts = [];

  if (any(tracking_options))
    for i=nframes:-1:2

      indx_interm = links{i}(:,1);
      indx_starts = [];

      nstarts = size(spots{i},1);
      indx_starts = setdiff([1:nstarts], indx_interm);
      nstarts = length(indx_starts);

      if (nstarts>0)
        starts(end+1:end+nstarts,:) = [spots{i}(indx_starts,:) indx_starts(:) ones(nstarts,1)*i];
      end

      nends = size(spots{i-1},1);
      indx_ends = setdiff([1:nends], links{i}(:,2));
      nends = length(indx_ends);
      if (nends>0)
        ends(end+1:end+nends,:) = [spots{i-1}(indx_ends,:) indx_ends(:) ones(nends,1)*i-1];
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEMS HERE !!!!!
      if (~isempty(spots{i}) & any(tracking_options(2:3)))
          
        join_dist = frame_weight(spots{i-1}(indx_ends,:), spots{i}(indx_interm,:));

        if (~isempty(prev_starts))
          split_dist = frame_weight(spots{i+1}(prev_starts, :), spots{i}(indx_interm, :));
          join_dist = [join_dist; split_dist];
        end

        indx_interm = indx_interm(any(join_dist < spot_max_movement, 1));
        ninterm = length(indx_interm);

        if (ninterm>0)
          interm(end+1:end+ninterm,:) = [spots{i}(indx_interm,:) indx_interm(:) ones(ninterm,1)*i];
        end
      end

      prev_starts = indx_starts;
    end
    
    nstarts = size(starts, 1);
    nends = size(ends, 1);
    ninterm = size(interm, 1);

    dist = Inf(nstarts + nends + 2*ninterm);

    if (tracking_options(1))
      mutual_dist = closing_weight(ends, starts);

      frame_indx = -bsxfun(@minus,ends(:,end),starts(:,end).');

      mutual_dist(frame_indx < 1 | frame_indx > max_frames | mutual_dist > spot_max_movement) = Inf;
    else
      mutual_dist = Inf(nends, nstarts);
    end

    if (tracking_options(2))
      [merge_dist, merge_weight, alt_weight] = joining_weight(ends, interm, spots, links);
      frame_indx = -bsxfun(@minus, ends(:,end), interm(:,end).');
      merge_weight = merge_dist .* merge_weight;
      merge_weight(merge_dist > spot_max_movement | frame_indx ~= 1) = Inf;
      alt_merge_weight = avg_movement * alt_weight;
    else
      merge_weight = Inf(nends, ninterm);
      alt_merge_weight = Inf(ninterm);
    end

    if (tracking_options(3))
      [split_dist, split_weight, alt_split_weight] = splitting_weight(interm, starts, spots, links);
      frame_indx = -bsxfun(@minus, interm(:,end), starts(:,end).');
      split_weight = split_dist .* split_weight;
      split_weight(split_dist > spot_max_movement | frame_indx ~= 1) = Inf;
      alt_split_weight = avg_movement * alt_split_weight;
    else
      split_weight = Inf(ninterm, nstarts);
      alt_split_weight = Inf(ninterm);
    end

    % Note that end-end merging and start-start splitting is not allowed by this
    % algorithm, which might make sense...

    dist(1:nends,1:nstarts) = mutual_dist;
    dist(1:nends,nstarts+1:nstarts+ninterm) = merge_weight;
    dist(nends+1:nends+ninterm,1:nstarts) = split_weight;

    trans_dist = dist(1:nends+ninterm, 1:nstarts+ninterm).';
    trans_dist(trans_dist < Inf) = min(trans_dist(:));

    alt_cost = prctile(dist(~isinf(dist)), 90);
    alt_dist = eye(max(nends, nstarts))*alt_cost;
    alt_dist(alt_dist==0) = Inf;

    dist(1:nends,nstarts+ninterm+1:end-ninterm) = alt_dist(1:nends,1:nends);
    dist(nends+ninterm+1:end-ninterm,1:nstarts) = alt_dist(1:nstarts,1:nstarts);

    dist(nends+1:nends+ninterm,end-ninterm+1:end) = alt_split_weight;
    dist(end-ninterm+1:end,nstarts+1:nstarts+ninterm) = alt_merge_weight;

    dist(nends+ninterm+1:end,nstarts+ninterm+1:end) = trans_dist;

    %disp('Let''s go !')

    %figure;
    %imagesc(dist)

    %[assign, cost] = munkres(dist);
    %[assign, cost] = assignmentoptimal(dist);
    [assign, cost] = lapjv_fast(dist);
    
  %  hold on;scatter(assign, [1:length(assign)])

    for i=1:nends+ninterm
      if (assign(i) <= nstarts)
        target = starts(assign(i), :);
      elseif (assign(i) < nstarts + ninterm)
        target = interm(assign(i) - nstarts, :);
      else
        continue;
      end

      if (i <= nends)
        %reference = [target(end-1) ends(i,end-1:end)];
        reference = ends(i,:);
      else
        %reference = [target(end-1) interm(i-nends,end-1:end)];
        reference = interm(i-nends,:);
      end

      if (opts.interpolate)
        ninterp = target(end) - reference(end);
        new_pts = bsxfun(@plus, bsxfun(@times, (reference(1:2) - target(1:2)) / ninterp, [1:ninterp-1].'), target(1:2));

        curr_pos = target(end-1);
        curr_indx = target(end);
        for j=1:ninterp-1
          curr_indx = target(end)-j;
          nprev = size(spots{curr_indx}, 1) + 1;
          spots{curr_indx} = [spots{curr_indx}; [new_pts(j,:) NaN(1,nprops)]];
          links{curr_indx+1} = [links{curr_indx+1}; [curr_pos nprev curr_indx]];
          curr_pos = nprev;
        end
        links{curr_indx} = [links{curr_indx}; [curr_pos, reference(end-1), curr_indx-1]];
      else
        links{target(end)} = [links{target(end)}; [target(end-1), reference(end-1:end)]];
      end
    end
  end

  if (opts.min_path_length > 0)
    min_length = opts.min_path_length;
    path_length = cell(nframes, 1);

    for i=1:nframes
      nimg = i;

      curr_links = links{nimg};
      path_length{nimg} = zeros(size(spots{nimg}, 1), 1);
      for j=1:size(curr_links, 1)
        path_length{nimg}(curr_links(j,1)) = path_length{curr_links(j,3)}(curr_links(j,2)) + nimg - curr_links(j,3);
      end
    end
    for i=nframes:-1:1
      nimg = i;

      curr_links = links{nimg};
      if (~isempty(curr_links))
        for j=1:size(curr_links, 1)
          path_length{curr_links(j,3)}(curr_links(j,2)) = path_length{nimg}(curr_links(j,1));
        end
      end

      long = (path_length{nimg} > min_length);
      good_indx = find(long);
      links{nimg} = curr_links(ismember(curr_links(:,1), good_indx), :);

      bad_indx = find(~long);
      spots{nimg}(bad_indx, :) = NaN;
      bad_prev = curr_links(ismember(curr_links(:,1), bad_indx), 2:3);
      if (~isempty(bad_prev))
        prev_frames = unique(bad_prev(:, 2)).';

        for p=prev_frames
          %interpolated = isnan(spots{nimg-1}(bad_prev, end));
          %spots{nimg-1}(bad_prev(interpolated), :) = NaN;
          curr_prev = bad_prev(bad_prev(:,2)==p,1);
          spots{p}(curr_prev, :) = NaN;
        end
      end
    end
  end

  if (~isempty(mymovie))
    if (isfield(mymovie, 'experiment'))
      if (opts.interpolate)
        for i=1:nframes
          mymovie.data.spots(i).cluster = links{i};
          mymovie.data.spots(i).carth = spots{i}(:, 1:2);
          mymovie.data.spots(i).properties = spots{i}(:, 3:end);
        end
      else
        for i=1:nframes
          mymovie.data.spots(i).cluster = links{i};
        end
      end
    else
      if (opts.interpolate)
        for i=1:nframes
          mymovie(i).cluster = links{i};
          if (~isempty(spots{i}))
            mymovie(i).carth = spots{i}(:, 1:2);
            mymovie(i).properties = spots{i}(:, 3:end);
          end
        end
      else
        for i=1:nframes
          mymovie(i).cluster = links{i};
        end
      end
    end
    links = mymovie;
  end

  return
end
