function mymovie = measure_flow(mymovie, opts)

  %opts = load_parameters(opts, 'track_flow');
  nframes = size_data(mymovie.data);

  if ((opts.recompute & opts.segment) | ~isfield(mymovie.data, 'spots') | isempty(mymovie.data.spots) | length(mymovie.data.spots)~=nframes | empty_struct(mymovie.data.spots, 'carth'))

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

  %{
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

  keyboard
  %}

  proj_dist = opts.spot_tracking.projection_dist / opts.pixel_size;
  pos_bin = opts.spot_tracking.projection_bin_size;
  frame_bin = opts.spot_tracking.projection_frames;
  if (isfield(opts.spot_tracking, 'projection_dist_min'))
    proj_min = opts.spot_tracking.projection_dist_min / opts.pixel_size;
  else
    proj_min = 0;
  end
  if (isfield(opts.spot_tracking, 'projection_args'))
    proj_args = opts.spot_tracking.projection_args;
  else
    proj_args = [];
  end
  if (strncmp(opts.spot_tracking.projection_type, 'gaussian', 8))
    if (isempty(proj_args))
      nullcline = 1.4;
    else
      nullcline = proj_args(1);
    end
    g_factor = 0.85;
  end

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

  nucleus_center = (isfield(mymovie.data, 'nuclei') & ~empty_struct(mymovie.data.nuclei, 'carth'));
  if (nucleus_center)
    for i=1:nframes
      if (~isempty(mymovie.data.nuclei(i).carth))
        prev_center = [mymovie.data.nuclei(i).carth mymovie.data.nuclei(i).properties i];
        break;
      end
    end
  end

  nucleus_width = 1.25;
  %nucleus_width = 3;

  all_dists = NaN(nframes, 2);

  init_shift = 1;
  %init_shift = 350

  prev_pts = mymovie.data.spots(init_shift).carth;
  for i=(init_shift+1):nframes
    nimg = i;
    %nimg = 450
    %nimg = i+350

    all_pts = mymovie.data.spots(nimg).carth;
    links = mymovie.data.spots(nimg).cluster;

    prev_links = links(links(:,end)==nimg-1, :);
    %good_links = (path_length{nimg}(prev_links(:,2)) >= path_thresh);
    %prev_links = prev_links(good_links, :);

    cortex = mymovie.data.cortex(nimg).carth;
    ncortex = size(cortex, 1);

    if (~isempty(prev_links) & ncortex > 0)

      pts = all_pts(prev_links(:,1), :);
      prev_pts = prev_pts(prev_links(:,2), :);

      npts = size(pts,1);

      egg_angle = mymovie.data.orientations(nimg);
      if (nucleus_center)
        if (~isempty(mymovie.data.nuclei(nimg).carth))
          nucleus_val = [mymovie.data.nuclei(nimg).carth mymovie.data.nuclei(nimg).properties nimg];
        else
          nucleus_val = prev_center;
        end

        ells = carth2elliptic([cortex; nucleus_val(1:2)], mymovie.data.centers(:, nimg), [1;1], 0, 'radial');
        ell_nucleus = ells(end,:);
        ells = ells(1:end-1, :);

        dangle = asin(nucleus_width*nucleus_val(3)/ell_nucleus(2));
  
        %ell_nucleus(1) = align_orientations(ell_nucleus(1), egg_angle);
        egg_angle = align_orientations(egg_angle, ell_nucleus(1));

        good_ells = (ells(:,1) >= ell_nucleus(1) - dangle & ells(:,1) <= ell_nucleus(1) + dangle);

        if (ell_nucleus(1) - dangle < 0)
          good_ells = good_ells | (ells(:,1)-2*pi >= ell_nucleus(1) - dangle & ells(:,1)-2*pi <= ell_nucleus(1) + dangle);
        end
        if (ell_nucleus(1) + dangle > 2*pi)
          good_ells = good_ells | (ells(:,1)+2*pi >= ell_nucleus(1) - dangle & ells(:,1)+2*pi <= ell_nucleus(1) + dangle);
        end

%        hold off;
%        plot(cortex(:,1), cortex(:,2), 'b');
%        hold on;
%        scatter(cortex(good_ells,1), cortex(good_ells,2), 'm');
%        rectangle
%        rectangle('Position', [nucleus_center(1:2)-nucleus_center(3) 2*nucleus_center(1,[3 3])], 'Curvature', [1 1], 'EdgeColor', 'k');
%        scatter(mymovie.data.centers(1,nimg), mymovie.data.centers(2,nimg), 'g');

%        drawnow

        dist = sqrt(mymean(sum(bsxfun(@minus, cortex(good_ells,:), nucleus_val(1:2)).^2, 2), 1));
        if (~isempty(dist))
          all_dists(i,1) = dist;
          dist = dist / nucleus_val(3);
          %rads & some function ....
          all_dists(i,2) = dist;

          %dangle = dangle*dist;
          % Sigmoid between 2 and 4 roughly, 3./(1+exp(-5*(x-2.5)))+1

            dangle = dangle*(4./(1+exp(-5*(dist-2.5)))+1);

          if (dist <= 2)
            egg_angle = ell_nucleus(1);
          elseif (dist <= 3)
            dist = dist - 2;
            egg_angle = ell_nucleus(1)*(1-dist) + egg_angle*dist;
          end

          good_ells = (ells(:,1) >= ell_nucleus(1) - dangle & ells(:,1) <= ell_nucleus(1) + dangle);

          if (ell_nucleus(1) - dangle < 0)
            good_ells = good_ells | (ells(:,1)-2*pi >= ell_nucleus(1) - dangle & ells(:,1)-2*pi <= ell_nucleus(1) + dangle);
          end
          if (ell_nucleus(1) + dangle > 2*pi)
            good_ells = good_ells | (ells(:,1)+2*pi >= ell_nucleus(1) - dangle & ells(:,1)+2*pi <= ell_nucleus(1) + dangle);
          end
        end

        prev_center = nucleus_val;
      end

      align_cortex = realign(cortex, [0;0], mymovie.data.centers(:, nimg), egg_angle);
      [pos, post_indxs, tot_dist] = carth2linear(align_cortex, true, opts);

%        hold off;
%        plot(align_cortex(:,1), align_cortex(:,2), 'b');
%        hold on;
%        scatter(align_cortex(good_ells,1), align_cortex(good_ells,2), 'm');
%        tmp_nucl = realign(nucleus_center(1:2), [0;0], mymovie.data.centers(:, nimg), egg_angle);
%        rectangle('Position', [tmp_nucl-nucleus_center(3) 2*nucleus_center(1,[3 3])], 'Curvature', [1 1], 'EdgeColor', 'k');
%        scatter(0, 0, 'g');
%
%        drawnow

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

          dist_thresh = nullcline / opts.pixel_size;

          weights = exp(-(dist.^2)/(2*((proj_dist/2)^2)));

          % Other means to compute the same thing, but less robust
          %{
          if (isfinite(dist_thresh))
            rescale = min(dist, [], 1);
            goods = (rescale < g_factor*dist_thresh);
            rescale = dist_thresh ./ (dist_thresh - rescale);
            weights(:,~goods) = 0;
            rescale(~goods) = 0;
          else
            rescale = ones(size(dist, 1), 1);
          end

          weights = bsxfun(@rdivide, weights, sum(weights, 2));
          weights = bsxfun(@times, weights, rescale);

          scaling = ones(size(dist, 1), 1);
          %}

          if (isfinite(dist_thresh))
            rescale = min(dist, [], 1);
            goods = (rescale < g_factor*dist_thresh);
            rescale(~goods) = 0;
            weights(:, ~goods) = 0;
            weights = bsxfun(@rdivide, weights, sum(weights, 2));
            rescale = bsxfun(@times, weights, rescale);
            rescale = sum(rescale, 2);
            goods = (rescale < g_factor*dist_thresh);
            scaling = dist_thresh ./ (dist_thresh - rescale);
            scaling(~goods) = NaN;
          else
            scaling = ones(size(dist, 1), 1);
          end

          weights = bsxfun(@rdivide, weights, sum(weights, 2));

          %{
          if (isfinite(dist_thresh))
            rescale = min(dist, [], 1);
            goods = (rescale < dist_thresh);
            rescale(~goods) = 0;
            weights(:, ~goods) = 0;
            weights = bsxfun(@rdivide, weights, sum(weights, 2));
            rescale = bsxfun(@times, weights, rescale);
            rescale = sum(rescale, 2);
            goods = (rescale < g_factor*dist_thresh);
            scaling = dist_thresh ./ (dist_thresh - rescale);
            scaling(~goods) = NaN;
          else
            scaling = ones(size(dist, 1), 1);
          end

          weights = bsxfun(@rdivide, weights, sum(weights, 2));
          %}

          speed_x = sum(bsxfun(@times, weights, speed(:,1).'), 2);
          speed_y = sum(bsxfun(@times, weights, speed(:,2).'), 2);

          movement = [speed_x, speed_y];
        case 'range'
          speed = pts(:,1:2) - prev_pts(:,1:2);
          dist = sqrt(bsxfun(@minus, cortex(:,1), pts(:,1).').^2 + ...
                      bsxfun(@minus, cortex(:,2), pts(:,2).').^2);

          close_ones = any(dist < proj_min, 1);
          weights = (dist <= proj_dist);
          weights(:, close_ones) = 0;
          good_ones = any(weights, 1);
          weights = bsxfun(@rdivide, weights, sum(weights, 2));
          %weights(isnan(weights)) = 0;

%          hold off;
%          plot(cortex(:,1), cortex(:,2), 'k');
%          hold on;
%          scatter(pts(:,1), pts(:,2), 'b');
%          scatter(pts(close_ones,1), pts(close_ones,2), 'r');
%          scatter(pts(good_ones,1), pts(good_ones,2), 'g');
%          drawnow

          scaling = ones(size(dist, 1), 1);
          %dist(isinf(dist)) = 0;
          %dist = dist .* weights;
          %dist = sum(dist, 2);
          %scaling = ones(size(dist));
          %scaling(dist < 0.75*dist_thresh) = dist_thresh ./ (dist_thresh - dist(dist < 0.75*dist_thresh));

          speed_x = sum(bsxfun(@times, weights, speed(:,1).'), 2);
          speed_y = sum(bsxfun(@times, weights, speed(:,2).'), 2);
%          speed_x = sum(bsxfun(@times, weights.*scaling, speed(:,1).'), 2);
%          speed_y = sum(bsxfun(@times, weights.*scaling, speed(:,2).'), 2);

          movement = [speed_x, speed_y];
        case 'nullcline'
          
          speed = pts(:,1:2) - prev_pts(:,1:2);

          speed_norm = sqrt(sum(speed.^2, 2));
          good_ones = (speed_norm < proj_min);

          dist = sqrt(bsxfun(@minus, cortex(:,1), pts(:,1).').^2 + ...
                      bsxfun(@minus, cortex(:,2), pts(:,2).').^2);

          dist(:, ~good_ones) = Inf;

%          hold off;
%          plot(cortex(:,1), cortex(:,2), 'k');
%          hold on;
%          scatter(pts(:,1), pts(:,2), 'b');
%          scatter(pts(good_ones,1), pts(good_ones,2), 'g');
%          drawnow


          scaling = min(dist, [], 2);

          movement = [perp(:,2), -perp(:,1)];
        otherwise
          error(['Projection type ''' opts.spot_tracking.projection_type ''' unknown.']);
      end

      if (nucleus_center)
   %     keyboard
        movement(good_ells(post_indxs), :) = NaN;
      end

      speed = dot([perp(:, 2), -perp(:, 1)], movement, 2) .* scaling;

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

      if (nimg == 290 && opts.verbosity == 3)
        img_size = [330 450]
        figure;imagesc(realign(double(load_data(mymovie.data, nimg)), img_size, mymovie.data.centers(:,nimg), mymovie.data.orientations(1,nimg)+pi));
        hold on;
        myplot(realign(cortex, img_size, mymovie.data.centers(:,nimg), mymovie.data.orientations(1,nimg)+pi));

        
        pts_tmp = (realign(pts(:,1:2), img_size, mymovie.data.centers(:,nimg), mymovie.data.orientations(1,nimg)+pi));
        scatter(pts_tmp(:,1), pts_tmp(:,2), 'k')

        prev_tmp = (realign(prev_pts(:,1:2), img_size, mymovie.data.centers(:,nimg), mymovie.data.orientations(1,nimg)+pi));
        scatter(prev_tmp(:,1), prev_tmp(:,2), 'y')
        speed_tmp = pts_tmp(:,1:2) - prev_tmp(:,1:2);

        plot([prev_tmp(:,1) pts_tmp(:,1)].', [prev_tmp(:,2) pts_tmp(:,2)].', 'w');

        quiver(pts_tmp(:,1), pts_tmp(:,2), speed_tmp(:,1), speed_tmp(:,2), 'k');

        curr_avg_flow = (realign([cortex movement], img_size, mymovie.data.centers(:,nimg), mymovie.data.orientations(1,nimg)+pi));
        quiver(curr_avg_flow(:,1), curr_avg_flow(:,2), curr_avg_flow(:,3), curr_avg_flow(:,4), 'r')

        curr_parall_flow = (realign([cortex perp(:, 2).*speed, -perp(:, 1).*speed], img_size, mymovie.data.centers(:,nimg), mymovie.data.orientations(1,nimg)+pi));
        quiver(curr_parall_flow(:,1), curr_parall_flow(:,2), curr_parall_flow(:,3), curr_parall_flow(:,4), 'g')
        keyboard;
      end
    else
      flow(nimg).speed = [];
      flow(nimg).position = [];
    end

    prev_pts = all_pts;
  end

  %keyboard

  mymovie.data.flow = flow;

  return;
end
