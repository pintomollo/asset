function mymovie = measure_flow(mymovie, opts)

  opts = load_parameters(opts, 'track_flow');
  nframes = size_data(mymovie.data);

  if (opts.recompute | ~isfield(mymovie.data, 'spots') | isempty(mymovie.data.spots) | length(mymovie.data.spots)~=nframes | empty_struct(mymovie.data.spots, 'carth'))

    mymovie = lidke_fit(mymovie, opts);
  end
  if (opts.recompute | empty_struct(mymovie.data.spots, 'cluster'))
    mymovie = track_spots(mymovie, opts);
  end

  if (~isfield(mymovie.data, 'flow'))
    flow = get_struct('flow', [1, nframes]);
  else  
    flow = mymovie.data.flow;
  end

  path_thresh = opts.spot_tracking.min_path_length;
  path_length = cell(nframes, 1);
  for i=1:nframes
    nimg = i;

    links = mymovie.data.spots(nimg).cluster;
    prev_links = links(links(:,end)==nimg-1, :);
    path_length{nimg} = zeros(size(mymovie.data.spots(nimg).carth, 1), 1);
    if (~isempty(prev_links))
      path_length{nimg}(prev_links(:,2)) = path_length{nimg-1}(prev_links(:,1)) + 1;
    end
  end
  for i=nframes:-1:1
    nimg = i;

    links = mymovie.data.spots(nimg).cluster;
    prev_links = links(links(:,end)==nimg-1, :);
    if (~isempty(prev_links))
      path_length{nimg-1}(prev_links(:,1)) = path_length{nimg}(prev_links(:,2));
    end
  end

  proj_dist = opts.spot_tracking.projection_dist / opts.pixel_size;
  pos_bin = opts.spot_tracking.projection_bin_size;
  frame_bin = opts.spot_tracking.projection_frames;

  edges = [0:pos_bin:65];
  edges = [-edges(end:-1:2) edges].';
  centers = edges(1:end-1) + pos_bin/2;
  edges(1) = -Inf;
  edges(end) = Inf;

  tmp_speed = NaN(size(centers));
  tmp_std = tmp_speed;

  all_speed = [];
  all_pos = [];
  all_frame = [];

  prev_pts = mymovie.data.spots(1).carth;
  for i=2:nframes
    nimg = i;

    all_pts = mymovie.data.spots(nimg).carth;
    links = mymovie.data.spots(nimg).cluster;

    prev_links = links(links(:,end)==nimg-1, :);
    good_links = (path_length{nimg}(prev_links(:,2)) >= path_thresh);
    prev_links = prev_links(good_links, :);

    cortex = mymovie.data.cortex(nimg).carth;
    ncortex = size(cortex, 1);

    if (~isempty(prev_links) & ncortex > 0)

      prev_pts = prev_pts(prev_links(:,1), :);
      pts = all_pts(prev_links(:,2), :);

      npts = size(pts,1);

      align_cortex = realign(cortex, [0;0], mymovie.data.centers(:, nimg), mymovie.data.orientations(nimg));
      [pos, post_indxs, tot_dist] = carth2linear(align_cortex, true, opts);
      cortex = cortex(post_indxs, :);
      [perp] = perpendicular_sampling(cortex, opts);
      pos = pos * tot_dist;

      switch (opts.spot_tracking.projection_type)
        case 'perp'
          vector_x = bsxfun(@minus, cortex(:,1), pts(:,1).');
          vector_y = bsxfun(@minus, cortex(:,2), pts(:,2).');

          dist = sqrt(vector_x.^2 + vector_y.^2);
          vector_x = vector_x ./ dist;
          vector_y = vector_y ./ dist;

          all_perp = repmat(perp, npts, 1);
          orient = acos(dot(all_perp, [vector_x(:) vector_y(:)], 2));
          orient = reshape(orient, ncortex, npts);

          goods = false(ncortex, npts);
          current = any(goods, 1);

          for p=5:5:50
            orient_thresh = prctile(orient(:, ~current), 2*p, 1);
            dist_thresh = prctile(dist(:, ~current), p, 1);

            goods(:, ~current) = bsxfun(@le, orient(:, ~current), orient_thresh) & bsxfun(@le, dist(:, ~current), dist_thresh);
            current = any(goods, 1);

            if (all(current))
              break;
            end
          end

          dist(~goods) = Inf;
          dist(dist>proj_dist) = Inf;

          [rel_dist, indxs] = min(dist, [], 1);
          good_pts = isfinite(rel_dist);
          indxs = indxs(good_pts);
          perp = perp(indxs, :);
          pos = pos(indxs);

          movement = pts(:,1:2) - prev_pts(:,1:2);
          movement = movement(good_pts,:);
        case 'gaussian'
          speed = pts(:,1:2) - prev_pts(:,1:2);
          dist = sqrt(bsxfun(@minus, cortex(:,1), pts(:,1).').^2 + ...
                      bsxfun(@minus, cortex(:,2), pts(:,2).').^2);
          dist(dist>3*proj_dist) = Inf;

          weights = exp(-(dist.^2)/(2*(proj_dist^2)));
          weights = bsxfun(@rdivide, weights, sum(weights, 2));

          speed_x = sum(bsxfun(@times, weights, speed(:,1).'), 2);
          speed_y = sum(bsxfun(@times, weights, speed(:,2).'), 2);

          movement = [speed_x, speed_y];
      end
      speed = dot([perp(:, 2), -perp(:, 1)], movement, 2);

      all_speed = [all_speed; speed*opts.pixel_size];
      all_pos = [all_pos; pos*opts.pixel_size];
      all_frame = [all_frame; nimg*ones(size(speed))];

      currents = (all_frame > nimg-frame_bin);
      all_speed = all_speed(currents);
      all_pos = all_pos(currents);
      all_frame = all_frame(currents);

      curr_speed = tmp_speed;
      curr_std = tmp_std;
      [counts, groups] = histc(all_pos, edges);
      [curr_speed(counts~=0), curr_std(counts~=0), groups] = mymean(all_speed, 1, groups);

      flow(nimg).speed = [curr_speed, curr_std];
      flow(nimg).position = centers;
    else
      flow(nimg).speed = [];
      flow(nimg).position = [];
    end

    prev_pts = all_pts;
  end

  mymovie.data.flow = flow;

  return;
end
