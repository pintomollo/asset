function ruffles = tracking(all_pts, opts)

  if (isstruct(all_pts))
    ruffles = all_pts;
    is_struct = true;
    nframes = length(ruffles);
  else
    nframes = max(all_pts(:,end));
    ruffles = get_struct('ruffles', [nframes 1]);
    is_struct = false;
  end

  max_thresh = 35;
  min_thresh = 5e-5;

  init_gap = 1.25;
  %gap_min = 100;

  window_size = 4;

  %if (opts.verbosity == 3)
  %  figure;hold on;
  %end

  prev_max = max_thresh;
  all_max = prev_max;

  list = [];
  pts = [];
  assign = [];
  npts = 0;
  avg_mov = [0 0];

  for i=1:nframes
    %if (i==76)
    %  keyboard
    %end
    %max_thresh = (mymovie.(type).axes_length(1,i) / 6);

    prev_pts = pts;
    prev_npts = npts;

    %pts = mymovie.(type).cortex(i).carth;
    %invs = find_ruffles(pts, mymovie.(type).centers(:,i),mymovie.(type).axes_length(:,i),mymovie.(type).orientations(1,i));

    %npts = size(invs, 1);

    %if (opts.verbosity == 3)
    %  ell_pos = carth2elliptic(mymovie.(type).cortex(i).carth, mymovie.(type).centers(:,i),mymovie.(type).axes_length(:,i),mymovie.(type).orientations(1,i));
    %  plot3(ell_pos(:,1),ell_pos(:,2),i*ones(size(ell_pos,1),1));
    %  hold on;
      %scatter3(ell_pos(invs(:,2),1),ell_pos(invs(:,2),2),i*ones(size(invs,1),1),'r')
    %end
  
    %pts = [pts(invs(:,2),:) pts(invs(:,1),:) pts(invs(:,3),:)];

    if (is_struct)
      pts = ruffles(i).carth;
    else
      pts = all_pts(all_pts(:,end) == i, 1:end-1);
    end
    npts = size(pts,1);
    
    %if (opts.verbosity == 3)
    %  scatter3(pts(:,1),pts(:,2),i*ones(npts,1),'r','LineWidth',2)
    %end

    if (prev_npts > 0 & npts > 0)

      dist = Inf(npts + prev_npts);
      
      mutual_dist = sqrt(bsxfun(@minus,prev_pts(:,1),pts(:,1).').^2 + bsxfun(@minus,prev_pts(:,2),pts(:,2).').^2);
      mutual_dist(mutual_dist > max_thresh) = Inf;

      if (i == 2)
        prev_max = max_thresh;
      end

      ends = eye(max(npts,prev_npts))*prev_max;
      ends(ends==0) = Inf;
      trans_dist = (mutual_dist.' < Inf) * min(mutual_dist(:));
      trans_dist(trans_dist == 0) = Inf;

      dist(1:prev_npts,1:npts) = mutual_dist;
      dist(1:prev_npts,npts+1:end) = ends(1:prev_npts,1:prev_npts);
      dist(prev_npts+1:end,1:npts) = ends(1:npts,1:npts);
      dist(prev_npts+1:end,npts+1:end) = trans_dist;

      %figure;implot(dist);

      [assign, cost] = munkres(dist);
      assign_dist = dist(sub2ind(size(dist),[1:prev_npts+npts],assign));
      good_indx = (assign(1:prev_npts) <= npts);
      avg_mov = avg_mov + [sum(assign_dist(good_indx)) sum(good_indx)];

      if (prev_max == max_thresh)
        prev_max = max(init_gap * assign_dist(assign_dist<prev_max));

        if (isempty(prev_max) | prev_max < min_thresh)
          prev_max = max_thresh;
        end

        all_max = prev_max;
      else
        prev_max = max([prev_max, init_gap * assign_dist(assign_dist<prev_max)]);
        all_max = max([all_max, assign_dist(assign_dist<prev_max)]);
      end
    end

    if (is_struct)
      ruffles(i) = store_detection(ruffles(i), prev_npts, assign);
    else
      ruffles(i) = store_detection(pts, prev_npts, assign);
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GAP CLOSING
  avg_mov = avg_mov(1) / avg_mov(2);

  starts = zeros(0, 3);
  ends = zeros(0, 4);
  interm = zeros(0, 5);

  for i=nframes:-1:2
    indxs = ruffles(i).cluster(:,1);
    new_starts = ruffles(i).carth(indxs == 0,:);
    nstarts = size(new_starts, 1);
    starts(end+1:end+nstarts,:) = [new_starts ones(nstarts,1)*i];

    nends = length(ruffles(i-1).cluster(:,1));
    indx_ends = setdiff([1:nends], indxs);
    nends = length(indx_ends);
    ends(end+1:end+nends,:) = [ruffles(i-1).carth(indx_ends,:) ones(nends,1)*i-1 ruffles(i-1).properties(indx_ends,:)];

    new_interm = find_intersect(ruffles, ends(end-nends+1:end,:), all_max, i);
    ninterm = size(new_interm, 1);
    interm(end+1:end+ninterm,:) = new_interm;
  end

  nstarts = size(starts, 1);
  nends = size(ends, 1);
  ninterm = size(interm, 1);

  dist = Inf(nstarts + nends + ninterm);

  mutual_dist = sqrt(bsxfun(@minus,ends(:,1),starts(:,1).').^2 + bsxfun(@minus,ends(:,2),starts(:,2).').^2);
  frame_indx = -bsxfun(@minus,ends(:,3),starts(:,3).');

  mutual_dist = mutual_dist ./ frame_indx;
  mutual_dist(frame_indx < 1 | frame_indx > window_size | mutual_dist > all_max) = Inf;

  alt_cost = prctile(mutual_dist(~isinf(mutual_dist)), 80);
  alt_dist = eye(max(nends, nstarts))*alt_cost;
  alt_dist(alt_dist==0) = Inf;

  if (isempty(mutual_dist))
    trans_dist = (mutual_dist.' < Inf);
  else
    trans_dist = (mutual_dist.' < Inf) * min(mutual_dist(:));
  end
  trans_dist(trans_dist == 0) = Inf;

  merge_dist = sqrt(bsxfun(@minus,ends(:,1),interm(:,1).').^2 + bsxfun(@minus,ends(:,2),interm(:,2).').^2);
  frame_indx = -bsxfun(@minus, ends(:,3), interm(:,3).');
  merge_dist(merge_dist > all_max | frame_indx ~= 1) = Inf;

  weight = repmat(interm(:,4).', nends, 1) ./ bsxfun(@plus,ends(:,4), interm(:,5).');
  weight(weight < 1) = weight(weight < 1).^(-2);
  merge_dist = merge_dist .* weight;

  alt_weight = diag((avg_mov*ones(ninterm,1)) .* interm(:,4) ./ interm(:,5));

  dist(1:nends,1:nstarts) = mutual_dist;
  dist(1:nends,nstarts+1:nstarts+ninterm) = merge_dist;
  dist(1:nends,nstarts+ninterm+1:end) = alt_dist(1:nends,1:nends);
  dist(nends+1:nends+ninterm,nstarts+1:nstarts+ninterm) = alt_weight(1:ninterm,1:ninterm);
  dist(nends+ninterm+1:end,1:nstarts) = alt_dist(1:nstarts,1:nstarts);
  dist(nends+ninterm+1:end,nstarts+ninterm+1:end) = trans_dist;

  [assign, cost] = munkres(dist);

  for i=1:nends
    if (assign(i) <= nstarts && assign(i) > 0)
      ruffles = close_gap(ruffles, ends(i,:), starts(assign(i),:));

      %if (opts.verbosity == 3)
      %  indx = ends(i, 3);
      %  ell_end = carth2elliptic(ends(i,1:2), mymovie.(type).centers(:,indx),mymovie.(type).axes_length(:,indx),mymovie.(type).orientations(1,indx));
      %  indx = starts(assign(i), 3);
      %  ell_start = carth2elliptic(starts(assign(i),1:2), mymovie.(type).centers(:,indx),mymovie.(type).axes_length(:,indx),mymovie.(type).orientations(1,indx));
      %  plot3([ell_end(1) ell_start(1)], [ell_end(2) ell_start(2)], [ends(i,3) starts(assign(i),3)], 'k');
      %end
    elseif (assign(i) < nstarts + ninterm && assign(i) > 0)
      ruffles = close_gap(ruffles, ends(i,:), interm(assign(i)-nstarts,:));

      %if (opts.verbosity == 3)
      %  indx = ends(i, 3);
      %  ell_end = carth2elliptic(ends(i,1:2), mymovie.(type).centers(:,indx),mymovie.(type).axes_length(:,indx),mymovie.(type).orientations(1,indx));
      %  indx = starts(assign(i), 3);
      %  ell_interm = carth2elliptic(interm(assign(i)-nstarts,1:2), mymovie.(type).centers(:,indx),mymovie.(type).axes_length(:,indx),mymovie.(type).orientations(1,indx));
      %  plot3([ell_end(1) ell_interm(1)], [ell_end(2) ell_interm(2)], [ends(i,3) interm(assign(i)-nstarts,3)], 'g');
      %end
    end
  end

  %mymovie.(type).ruffles = ruffles;

  %if (opts.verbosity == 3)
  %  view(0,0);
  %end

  return
end

function merges = find_intersect(ruffles, ends, thresh, frame)

  merges = zeros(0, 4);

  valids = (ruffles(frame).cluster(:,1) ~= 0);
  indxs = find(valids);

  interms = ruffles(frame).carth(valids,:);

  if (isempty(interms))
    return;
  end

  mutual_dist = sqrt(bsxfun(@minus,ends(:,1),interms(:,1).').^2 + bsxfun(@minus,ends(:,2),interms(:,2).').^2);

  cands = indxs(any(mutual_dist < thresh,1));
  ncands = length(cands);
  prev_indx = ruffles(frame).cluster(cands,1);

  merges = [ruffles(frame).carth(cands,:) ones(ncands,1)*frame ruffles(frame).properties(cands,:) ruffles(frame-1).properties(prev_indx,:)];

  return;
end

function ruffles = close_gap(ruffles, first, last)

  %keyboard

  first_indx = first(3);
  last_indx = last(3);

  indx = find(all(ruffles(first_indx).carth(:,1) == first(1,1) & ruffles(first_indx).carth(:,2) == first(1,2), 2));
  indx2 = find(all(ruffles(last_indx).carth(:,1) == last(1,1) & ruffles(last_indx).carth(:,2) == last(1,2), 2));

  dist = last_indx - first_indx;

  if (dist > 1)

    dmov = diff([first(1,1:2); last(1,1:2)]);
    interm_pos = first(ones(1,dist-1),1:2) +  ([1:dist-1] / dist).' * dmov;

    for i=1:dist-1
      frame_indx = i+first_indx;
      ruffles(frame_indx).carth(end+1,:) = interm_pos(i,:);
      ruffles(frame_indx).cluster(end+1,1) = indx;
      %ruffles(frame_indx).bounds(end+1,:) = 0;
      ruffles(frame_indx).properties(end+1,:) = 0;

      indx = length(ruffles(frame_indx).cluster(:,1));
    end

    ruffles(last_indx).cluster(indx2,1) = indx;
  else
    if (ruffles(last_indx).cluster(indx2,1) == 0)
      ruffles(last_indx).cluster(indx2,1) = indx;
    else
      ruffles(last_indx).cluster(indx2,2) = indx;
    end
  end

  return;
end

function ruffles = store_detection(pts, nprev, assign)

  if (isstruct(pts))
    ruffles = pts;
    npts = size(ruffles.carth, 1);
  else

    npts = size(pts,1);

    %ruffles = get_struct('ruffles', 1); 
    %ruffles.carth = pts(:,1:2);
    %ruffles.bounds = pts(:,3:end);

    ruffles = get_struct('ruffles');
    ruffles.carth = pts(:,1:2);
    ruffles.properties = pts(:,3:end);
  end

  if (nprev > 0)
    [junk tmp] = sort(assign);
      
    ruffles.cluster = [tmp(1:npts).' zeros(npts, 1)];
    ruffles.cluster(ruffles.cluster > nprev,1) = 0;
  else
    ruffles.cluster = zeros(npts, 2);
  end

  %figure;
  %hold on;
  %patch(pts(:,[1:2:end]).', pts(:,[2:2:end]).','r')
  %keyboard

  return
end
