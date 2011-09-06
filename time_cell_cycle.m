function mymovie = time_cell_cycle(mymovie, opts)

  type = opts.segmentation_type;
  if (strncmp(type, 'markers', 7))
    [nframes imgsize] = size_data(mymovie.cortex);
  else
    [nframes imgsize] = size_data(mymovie.(type));
  end

  if (~opts.recompute && isfield(mymovie.(type), 'timing') && ~isempty(mymovie.(type).timing) && ~isnan(mymovie.(type).timing.cytokinesis))
    
    return;
  end

  if (~isfield(mymovie.(type), 'timing'))
    mymovie.(type).timing = get_struct('timing');
  end

  if (~isfield(mymovie.(type), 'ruffles')|isempty(mymovie.(type).ruffles)|~isfield(mymovie.(type).ruffles, 'paths')|isempty(mymovie.(type).ruffles(1).paths)|isempty(mymovie.(type).ruffles(1).cluster))
    mymovie = track_ruffles(mymovie, opts);
  end

  if (~isfield(mymovie.dic, 'nuclei') | isempty(mymovie.dic.nuclei))
    mymovie = detect_dic_nuclei(mymovie, opts);
  end
 
  nbins = 16;

  cortex = cell(1, nframes);
  pts = zeros(0, 3);

  for i=1:nframes
    carth = mymovie.(type).ruffles(i).carth;
    ell = carth2elliptic(carth, mymovie.(type).centers(:,i),mymovie.(type).axes_length(:,i),mymovie.(type).orientations(1,i));
    ell = [ell i*ones(size(ell, 1), 1)];
    pts = [pts; ell];
    cortex{i} = carth2elliptic(mymovie.(type).cortex(i).carth, mymovie.(type).centers(:,i),mymovie.(type).axes_length(:,i),mymovie.(type).orientations(1,i));
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
        traces{groups(j, 1)} = [traces{groups(j, 1)}; [ell_pos(groups(j, 2),:) mymovie.(type).ruffles(i).properties(groups(j,2), :)]];
        groups(j, 2) = mymovie.(type).ruffles(i).cluster(groups(j,2), 1);
      end
      groups = groups(groups(:,2) ~= 0, :);

      frames(i) = true;
    end

    ell_bin = indx(pts(:,3)==i);
    valids = true(size(ell_pos(:,1)));
    %valids = (ell_bin == cytok_pos(1)|ell_bin == cytok_pos(2));
    if (~isempty(done_pts))
      valids(done_pts) = false;
    end

    new_indx = (valids & mymovie.(type).ruffles(i).cluster(:,1) ~= 0);
    new_pts = [ell_pos(new_indx, :) mymovie.(type).ruffles(i).properties(new_indx, :)];
    new_cluster = mymovie.(type).ruffles(i).cluster(new_indx, 1);
    for j=1:length(new_cluster)
      traces{end+1} = new_pts(j, :);
      groups = [groups; [length(traces) new_cluster(j)]];
    end

    if (~frames(i) & ~isempty(new_pts))
      frames(i) = true;
    end
  end

  traces{end+1} = zeros(0, 5);

  params = [0.7 0.2 0.8];

  [left_cytok, right_cytok] = find_pairs(frames, traces, params);
  cytok_frame = min([traces{left_cytok}(:, 3); traces{right_cytok}(:, 3)]);
  cytok_pos = median([traces{left_cytok}(:, 1); 2*pi - traces{right_cytok}(:, 1)]);

  if (cytok_frame == 1)
    mymovie.metadata.timing.cytokinesis = NaN;
    cytok_frame = nframes;
  else
    mymovie.metadata.timing.cytokinesis = cytok_frame;
  end

  params = [0.5 0.7 0.5];
  pc_frames = frames;
  pc_frames(find(pc_frames, 1, 'first'):find(pc_frames(1:cytok_frame-1), 1, 'last')) = true;
  [left_pc, right_pc] = find_pairs(pc_frames(1:cytok_frame-1), traces, params, cytok_pos);

  max_frame = max(traces{left_pc}(:,3));
  [lambda_left, alpha_left] = estimate_gamma(max_frame - traces{left_pc}(:,3) + 1, 1 - traces{left_pc}(:, 2));
  [lambda_left, alpha_left] = estimate_gamma(max_frame - traces{left_pc}(:,3) + 1, traces{left_pc}(:, 5));

  pc_left = round(max_frame - ((alpha_left - 1) / lambda_left) + 1);
  max_frame = max(traces{right_pc}(:,3));
  [lambda_right, alpha_right] = estimate_gamma(max_frame - traces{right_pc}(:,3) + 1, 1 - traces{right_pc}(:, 2));
  [lambda_right, alpha_right] = estimate_gamma(max_frame - traces{right_pc}(:,3) + 1, traces{right_pc}(:, 5));
  pc_right = round(max_frame - ((alpha_right - 1) / lambda_right) + 1);

  %keyboard

  if (false)
  figure;hold on;
  for i=1:nframes
    plot3(cortex{i}(:,1),cortex{i}(:,2),i*ones(size(cortex{i}(:,2))), 'b');
  end
  for i=1:length(traces);
    if (i == left_pc|i == right_pc)
      colors = 'r';
    elseif (i == left_cytok|i == right_cytok)
      colors = 'g';
    else
      colors = 'k';
    end
    plot3(traces{i}(:,1), traces{i}(:,2), traces{i}(:,3), ['-o' colors]);
  end
  end

  pc_frame = ceil(mymean([pc_left; pc_right]));

  mymovie.metadata.timing.pseudocleavage = pc_frame;

  traces = {};
  prev_pts = [];
  groups = [];

  for i=nframes:-1:1
  
    pts = [mymovie.dic.nuclei(i).carth mymovie.dic.nuclei(i).properties];

    if (isempty(pts))
      continue;
    end

    done_pts = [];
    if (~isempty(groups))
      done_pts = groups(:,2);
      for j=1:size(groups, 1)
        traces{groups(j, 1)} = [traces{groups(j, 1)}; [pts(groups(j, 2),:) i]];
        groups(j, 2) = mymovie.dic.nuclei(i).cluster(groups(j,2), 1);
      end
      groups = groups(groups(:,2) ~= 0, :);
    end

    valids = true(size(pts(:,1)));
    if (~isempty(done_pts))
      valids(done_pts) = false;
    end

    new_indx = (valids & mymovie.dic.nuclei(i).cluster(:,1) ~= 0);
    new_pts = pts(new_indx, :);
    new_cluster = mymovie.dic.nuclei(i).cluster(new_indx, 1);
    for j=1:length(new_cluster)
      traces{end+1} = [new_pts(j, :) i];
      groups = [groups; [length(traces) new_cluster(j)]];
    end
  end

  lengths = zeros(length(traces), 1);
  for i=1:length(traces)
    max_frame = max(traces{i}(:,end));
    min_frame = min(traces{i}(:,end));

    if (max_frame >= pc_frame & min_frame <= cytok_frame)
      lengths(i) = max_frame - min_frame + 1;
    end
  end

  traces = traces(lengths > 0);
  lengths = lengths(lengths > 0);

  traces{end+1} = zeros(0, 4);
  traces{end+1} = zeros(0, 4);
  lengths(end+1) = 0;
  lengths(end+1) = 0;

  [~, indxs] = sort(lengths);
  bests = indxs([end end-1]);
  traces = traces(bests);
  nucl1 = traces{1};
  nucl2 = traces{2};
  
  thresh = 10;
  pnm_frame = NaN;
  for i=pc_frame:cytok_frame
    pt1 = nucl1(nucl1(:,end) == i, :);
    pt2 = nucl2(nucl2(:,end) == i, :);

    if (isempty(pt1) | isempty(pt2))
      continue;
    end

    dist = sqrt(sum((pt1(1:2) - pt2(1:2)).^2));
    radii = sqrt(pt1(3)/pi) + sqrt(pt2(3)/pi);

    if (dist < radii + thresh)
      pnm_frame = i;
      break;
    end
  end

  if (pnm_frame == pc_frame)
    pnm_frame = NaN;
  end

  mymovie.metadata.timing.pronuclear_meeting = pnm_frame;

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

  %%figure;
  %hold on;

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

    %plot(traces{i}(:,3), traces{i}(:,1));
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

  %keyboard

  return;
end
