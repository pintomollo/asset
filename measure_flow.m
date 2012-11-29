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

  prev_pts = mymovie.data.spots(1).carth;
  for i=2:nframes
    nimg = i;

    all_pts = mymovie.data.spots(nimg).carth;
    links = mymovie.data.spots(nimg).cluster;

    prev_links = links(links(:,end)==nimg-1, :);
    cortex = mymovie.data.cortex(nimg).carth;
    ncortex = size(cortex, 1);

    if (~isempty(prev_links) & ncortex > 0)

      prev_pts = prev_pts(prev_links(:,1), :);
      pts = all_pts(prev_links(:,2), :);

      npts = size(pts,1);

      [perp, junk, linear] = perpendicular_sampling(cortex, opts);

      vector_x = bsxfun(@minus, cortex(:,1), pts(:,1).');
      vector_y = bsxfun(@minus, cortex(:,2), pts(:,2).');

      dist = sqrt(vector_x.^2 + vector_y.^2);
      vector_x = vector_x ./ dist;
      vector_y = vector_y ./ dist;

      all_perp = repmat(perp, npts, 1);
      orient = acos(dot(all_perp, [vector_x(:) vector_y(:)], 2));
      %orient = min(abs(orient), (pi - orient));
      orient = reshape(orient, ncortex, npts);

      %angle_thresh = pi/16;
      prcnt_thresh = 1;
      goods = false(ncortex, npts);
      current = any(goods, 1);

      for p=5:5:50
      %while(~all(current))
        orient_thresh = prctile(orient(:, ~current), 2*p, 1);
        dist_thresh = prctile(dist(:, ~current), p, 1);

        %goods(:, ~current) = (orient(:, ~current) < angle_thresh);
        goods(:, ~current) = bsxfun(@le, orient(:, ~current), orient_thresh) & bsxfun(@le, dist(:, ~current), dist_thresh);
        current = any(goods, 1);

        if (all(current))
          %display(p)
          break;
        end
        %angle_thresh = angle_thresh + angle_thresh;
      end

      %img = imnorm(double(load_data(mymovie.data, nimg)));
      %hold off;
      %imshow(img);
      %hold on;
      %scatter(pts(:,1), pts(:,2), 'r');
      %quiver(pts(:,1), pts(:,2), pts(:,1) - prev_pts(:,1), pts(:,2) - prev_pts(:,2), 0, 'r');
      %print('-dpng', '-r150', ['PNG/Flow/tracking_' num2str(nimg) '.png']);

      %cortex = realign(cortex, [0;0], mymovie.data.centers(:, nimg), mymovie.data.orientations(nimg));
      %hold off;
      %plot(cortex(:,1), cortex(:,2), 'k')
      %hold on;
      %quiver(pts(:,1), pts(:,2), pts(:,1) - prev_pts(:,1), pts(:,2) - prev_pts(:,2), 'b');
      %drawnow
      %pause(0.25);

      dist(~goods) = Inf;
      [rel_dist, indxs] = min(dist, [], 1);
      rel_perp = perp(indxs, :);
      rel_pos = linear(indxs);

      proj_speed = dot([rel_perp(:, 2), -rel_perp(:, 1)], pts(:,1:2) - prev_pts(:,1:2), 2);

      flow(nimg).speed = proj_speed;
      flow(nimg).distance = rel_dist;
      flow(nimg).position = rel_pos;
      flow(nimg).index = indxs;
      %keyboard
    else
      flow(nimg).speed = [];
      flow(nimg).distance = [];
      flow(nimg).position = [];
      flow(nimg).index = [];
    end

    prev_pts = all_pts;
  end

  mymovie.data.flow = flow;

  return;
end
