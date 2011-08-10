function [pc_frame, cytok_frame] = time_cell_cycle(mymovie, opts)
  
  [nframes] = size_data(mymovie.dic);
  nbins = 16;

  pts = zeros(0, 3);

  for i=1:nframes
    carth = mymovie.markers.ruffles(i).carth;
    ell = carth2elliptic(carth, mymovie.markers.centers(:,i),mymovie.markers.axes_length(:,i),mymovie.markers.orientations(1,i));
    ell = [ell i*ones(size(ell, 1), 1)];
    pts = [pts; ell];
  end

  bins = [0:2*pi/nbins:2*pi + 1e-5];
  [counts, indx] = histc(pts(:,1), bins);

  myhist = counts(1:nbins/2) + counts(end-1:-1:end-nbins/2);

  size_front = nbins / 8;
  [~, cytok_pos] = max(myhist(size_front+1:end-size_front));
  cytok_pos = cytok_pos + size_front;
  cytok_pos = [cytok_pos nbins - cytok_pos + 1];


  traces = {};
  prev_pts = [];
  groups = [];
  frames = false(1, nframes);

  for i=nframes:-1:1
  
    ell_pos = pts(pts(:,3)==i, :);

    done_pts = [];
    if (~isempty(groups))
      done_pts = groups(:,2);
      for j=1:size(groups, 1)
        traces{groups(j, 1)} = [traces{groups(j, 1)}; ell_pos(groups(j, 2),:)];
        groups(j, 2) = mymovie.markers.ruffles(i).cluster(groups(j,2), 1);
      end
      groups = groups(groups(:,2) ~= 0, :);

      frames(i) = true;
    end

    ell_bin = indx(pts(:,3)==i);
    valids = (ell_bin == cytok_pos(1)|ell_bin == cytok_pos(2));
    if (~isempty(done_pts))
      valids(done_pts) = false;
    end

    new_indx = (valids & mymovie.markers.ruffles(i).cluster(:,1) ~= 0);
    new_pts = ell_pos(new_indx, :);
    new_cluster = mymovie.markers.ruffles(i).cluster(new_indx, 1);
    for j=1:length(new_cluster)
      traces{end+1} = new_pts(j, :);
      groups = [groups; [length(traces) new_cluster(j)]];
    end

    if (~frames(i) & ~isempty(new_pts))
      frames(i) = true;
    end
  end

  traces{end+1} = zeros(0, 3);

  params = [0.7 0.2 0.8];

  [left_cytok, right_cytok] = find_pairs(frames, traces, params);
  cytok_frame = min([traces{left_cytok}(:, 3); traces{right_cytok}(:, 3)]);
  cytok_pos = median([traces{left_cytok}(:, 1); 2*pi - traces{right_cytok}(:, 1)]);

  params = [0.5 0.7 0.5];
  pc_frames = frames;
  pc_frames(find(pc_frames, 1, 'first'):find(pc_frames(1:cytok_frame-1), 1, 'last')) = true;
  [left_pc, right_pc] = find_pairs(pc_frames(1:cytok_frame-1), traces, params, cytok_pos);


  %pc_position(traces{left_pc}(:, [3 2]), traces{right_pc}(:, [3 2]));

   %[both_x, ~, groups] = unique([traces{left_pc}(:,3); traces{right_pc}(:,3)]);
   % y = mymean([traces{left_pc}(:,2); traces{right_pc}(:,2)], [], groups);

   % x = both_x;
   % max_frame = max(x);
   % x = max(x) - x + 1;
   % y = 1 - y;

   % [lambda, alpha] = estimate_gamma(x, y);
   [lambda_left, alpha_left] = estimate_gamma(max(traces{left_pc}(:,3)) - traces{left_pc}(:,3) + 1, 1 - traces{left_pc}(:, 2));
    pc_left = round(max_frame - ((alpha - 1) / lambda) + 1);
   [lambda_right, alpha_right] = estimate_gamma(max(traces{right_pc}(:,3)) - traces{right_pc}(:,3) + 1, 1 - traces{right_pc}(:, 2));
    pc_right = round(max_frame - ((alpha - 1) / lambda) + 1);

    pc_frame = mymean([pc_left, pc_right]);


  %figure;scatter3(pts(:,1), pts(:,2), pts(:,3)); hold on;
  %plot3(traces{left_cytok}(:,1), traces{left_cytok}(:,2), traces{left_cytok}(:,3), 'r');
  %plot3(traces{right_cytok}(:,1), traces{right_cytok}(:,2), traces{right_cytok}(:,3), 'r');
  %plot3(traces{left_pc}(:,1), traces{left_pc}(:,2), traces{left_pc}(:,3), 'g');
  %plot3(traces{right_pc}(:,1), traces{right_pc}(:,2), traces{right_pc}(:,3), 'g');

  return;
end

function [left_indx, right_indx] = find_pairs(frames, traces, params, target)

  if (nargin < 4)
    target = pi/2;
  end

  alpha = params(1);
  beta = params(2);
  gamma = params(3);

  frame = 0;
  pos = 0;

  nframes = length(frames);

  max_cytok = find(frames, 1, 'last');
  if (isempty(max_cytok))
    max_cytok = 1;
  end
  min_cytok = find(~frames(1:max_cytok), 1, 'last');
  if (isempty(min_cytok))
    min_cytok = 1;
  end

  left = [];
  right = [];
  for i=1:length(traces)-1
    if (max(traces{i}(:,3)) < min_cytok | min(traces{i}(:,3)) > max_cytok)
      continue;
    end
    
    tmp_pos = median(traces{i}(:,1));
    if (tmp_pos < pi)
      left = [left; [tmp_pos, size(traces{i}, 1), i, max(traces{i}(:,3)), min(traces{i}(:, 3))]];
    else
      right = [right; [tmp_pos, size(traces{i}, 1), i, max(traces{i}(:,3)), min(traces{i}(:, 3))]];
    end
  end

  dist = [];
  left_indx = 0;
  right_indx = 0;
  norm_factor = (max_cytok - min_cytok + 1);

  if (numel(left) == 0)
    if (numel(right) ~= 0)
      dist = alpha * abs(target - (2*pi - right(:,1))) / (pi/8) + ...
             (1-alpha)*(beta * (1 - (right(:,2) / norm_factor)) + ...
              (1-beta) * (nframes - right(:, 4)) / norm_factor);
      
      [~, right_indx] = min(dist);
    end
  elseif (numel(right) == 0)
    dist = alpha * abs(left(:, 1) - target) / (pi/8) + ...
           (1-alpha)*(beta * (1 - (left(:,2) / norm_factor)) + ....
            (1-beta) * (nframes - left(:, 4)) / norm_factor);

    [~, left_indx] = min(dist);
  else
    dist = Inf(size(left, 1), size(right, 1));
    for i=1:size(left, 1)
      for j=1:size(right, 1)
        if (left(i, 4) < right(j, 5) | left(i, 5) > right(j, 4))
          continue;
        end

        dist(i,j) = alpha * (beta* ((1 - (left(i,2) / norm_factor)) + ....
                    (1 - (right(j,2) / norm_factor))) / 2 + ...
                    (1-beta) * ((nframes - left(i, 4)) / norm_factor + ...
                    (nframes - right(j, 4)) / norm_factor) / 2) + ...
                    (1-alpha) * (gamma* abs(left(i, 1) - (2*pi - right(j,1))) / (pi/4) + ...
                    (1-gamma) * (1 - (min(right(j,4), left(i, 4)) - max(right(j, 5), left(i,5))) / norm_factor));
      end 
    end
    tmp = min(dist(:));
    [left_indx, right_indx] = find(dist == tmp);
  end
        %dist(i,j) = alpha * abs(left(i, 1) - (2*pi - right(j,1))) / (pi/4) + ...
        %            (1-alpha)*(beta * ((1 - (left(i,2) / norm_factor)) + ....
        %            (1 - (right(j,2) / norm_factor))) / 2 + ...
        %            (1-beta) * ((nframes - left(i, 4)) / norm_factor + ...
        %            (nframes - right(j, 4)) / norm_factor) / 2);

        %dist(i,j) = alpha * (beta* abs(left(i, 1) - (2*pi - right(j,1))) / (pi/4) + ...
        %            (1-beta)* ((1 - (left(i,2) / norm_factor)) + ....
        %            (1 - (right(j,2) / norm_factor)))/2) + ....
        %            (1-alpha)*(gamma * (abs(right(j,2) - left(i,2)) / norm_factor) + ....
        %            (1-gamma) * (1 - (min(right(j,4), left(i, 4)) - max(right(j, 5), left(i,5))) / norm_factor));
  %keyboard

  if (left_indx ~= 0)
    left_indx = left(left_indx, 3);
  else
    left_indx = length(traces);
  end
  if (right_indx ~= 0)
    right_indx = right(right_indx, 3);
  else
    right_indx = length(traces);
  end

  return;
end
