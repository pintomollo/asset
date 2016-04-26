%function [flow, flow_std, pos] = display_flow(mymovie, opts, args)
function [flow, flow_std, pos] = display_flow(mymovie, opts, do_piv)

  if (nargin < 3)
    do_piv = false;
  end

  if (do_piv)
    mymovie.data.flow = mymovie.data.piv;
  end

%  thresh = opts.quantification.pole_threshold;
  nframes = size_data(mymovie.data);
  
%  if (nargin < 3)
%    dist_thresh = 1;
%    pos_bin = 5;
%    frame_bin = 10;
%  else
%    dist_thresh = args(1);
%    pos_bin = args(2);
%    frame_bin = args(3);
%  end

%  edges = [0:pos_bin:50];
%  edges = [-edges(end:-1:2) edges].';
%  centers = edges(1:end-1) + pos_bin/2;
%  edges(1) = -Inf;
%  edges(end) = Inf;

  flow = [];
  flow_std = [];
  pos = [];

  if (empty_struct(mymovie.data.flow, 'position'))
    return;
  end

  centers = [];
  for i=1:nframes
    centers = mymovie.data.flow(i).position;
    if (~isempty(centers))
      break;
    end
  end

%  all_speed = NaN(size(centers));
%  all_dist = zeros(size(centers));
%  all_pos = centers;
%  all_frame = Inf(size(centers));

  if (nargout > 0)
    flow = NaN(length(centers), nframes);
    flow_std = flow;
  end

  if (~isfield(opts.spot_tracking, 'scaling_parameter'))
    opts = merge_structures(opts, get_struct('ASSET'));
  end
  if (do_piv)
    scaling = 1;
  else
    scaling = opts.spot_tracking.scaling_parameter;
  end

  for i=1:nframes
    nimg = i;

    %speed = mymovie.data.flow(nimg).speed * opts.pixel_size;
    %dist = mymovie.data.flow(nimg).distance * opts.pixel_size;
    %pos = mymovie.data.flow(nimg).position;
    speed = mymovie.data.flow(nimg).speed * scaling;
    %indxs = mymovie.data.flow(nimg).index;

    %cortex = mymovie.data.cortex(nimg).carth;
    
    %if (isempty(cortex) | isempty(speed))
    if (isempty(speed))
      continue;
    end

    curr_speed = speed(:,1);
    curr_std = speed(:,2);

    %{
    cortex = realign(cortex, [0;0], mymovie.data.centers(:, nimg), mymovie.data.orientations(nimg));
    cortex = cortex * opts.pixel_size;

    [pts, indexes, dists] = realign_poles(cortex, thresh);

    [junk, revert] = sort(indexes);
    pts = pts(revert);
    pts = pts(indxs);
    %goods = ismember(indexes, indxs);
    %if (any(goods))

    goods = dist < dist_thresh;
    pts = pts(goods);
    dist = dist(goods).';
    speed = speed(goods);

    all_speed = [all_speed; speed];
    all_dist = [all_dist; dist];
    all_pos = [all_pos; pts];
    all_frame = [all_frame; nimg*ones(size(speed))];

    currents = (all_frame > nimg-frame_bin);
    all_speed = all_speed(currents);
    all_dist = all_dist(currents);
    all_pos = all_pos(currents);
    all_frame = all_frame(currents);

    [counts, groups] = histc(all_pos, edges);
    [curr_speed, curr_std, groups] = mymean(all_speed, 1, groups);

    %keyboard

    %curr_speed = curr_speed(groups ~= 0);
    %curr_std = curr_std(groups ~= 0);

    %figure;quiver(pos, dist, speed, zeros(size(speed)));

    %  tmp_indxs = indexes(goods);
    %  pts = pts(goods);

    %  [vals, revert] = sort(tmp_indxs);
    %  speed = speed(revert);
    %  dist = dist(revert);
%    if (mod(nimg, frame_bin) == 0)
    %}

    if (nargout == 0)
      subplot(1,2,1)
      hold off;errorbar(centers, curr_speed, curr_std);
      subplot(1,2,2)
      hold off;
      imshow(imnorm(double(load_data(mymovie.data, nimg))));
      %hold off;quiver(pts(goods), dist(goods), speed(goods), zeros(sum(goods), 1));
      drawnow
      pause(0.25);
    else
      flow(:, nimg) = curr_speed;
      flow_std(:, nimg) = curr_std;
    end
%    end
  end

  if (nargout == 2)
    flow_std = centers;
  else
    pos = centers;
  end

  return;
end

function [pts, indexes, dists] = realign_poles(cortex, thresh)

  [pts, total_dist] = carth2linear(cortex);
  
  goods = (cortex(:,1) > 0);
  max_pos = max(cortex(goods, 1));
  goods = (cortex(:,1) > (1-thresh)*max_pos);

  max_r = min(abs(cortex(goods,2)));
  post_indx = find(goods & abs(cortex(:,2)) == max_r, 1);

  goods = (cortex(:,1) < 0);
  max_pos = min(cortex(goods, 1));
  goods = (cortex(:,1) < (1-thresh)*max_pos);

  max_r = min(abs(cortex(goods,2)));
  ant_indx = find(goods & abs(cortex(:,2)) == max_r, 1);

  pts = pts * total_dist;
  pts = pts - pts(post_indx);
  if (ant_indx > post_indx)
    pts = [pts(ant_indx:end, :) - total_dist; pts(1:ant_indx-1, :)];
  else
    pts = [pts(ant_indx:end, :); pts(1:ant_indx-1, :) + total_dist];
  end

  dists = [pts(1), total_dist + pts(1)];
  indexes = [ant_indx:length(pts) 1:ant_indx-1];

  return;
end
