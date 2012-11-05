function mymovie = find_ruffles(mymovie, opts)
%pts, center, axes_length, orient, hard_thresh)

  if (nargin < 2)
    opts = get_struct('ASSET', 1);
  end

  type = opts.segmentation_type;
  switch (type)
    case 'all'
      [nframes imgsize ] = size_data(mymovie.dic);
      if (isfield(mymovie, 'markers' & ~isempty(mymovie.markers)))
        type = 'markers';
      else
        type = 'dic';
      end
    case 'dic'
      [nframes imgsize ] = size_data(mymovie.dic);
    case 'markers'
      if (isfield(mymovie, 'eggshell') & ~isempty(mymovie.eggshell))
        [nframes imgsize ] = size_data(mymovie.eggshell);
      else
        [nframes imgsize ] = size_data(mymovie.cortex);
      end
    case 'data'
      [nframes imgsize ] = size_data(mymovie.data);
    otherwise
      error 'None of the expected field are present in ''mymovie''';
  end

  if (nargin < 2)
    opts.segmentation_parameters = set_image_size(opts.segmentation_parameters, imgsize);
  end

  if (~opts.recompute && isfield(mymovie.(type), 'ruffles') && ~isempty(mymovie.(type).ruffles) && ~any(isnan(mymovie.(type).ruffles(1).carth(:))))
    return;
  end

  if (~isfield(mymovie.(type), 'cortex') || isempty(mymovie.(type).cortex))
    mymovie = segment_movie(mymovie, opts);
  end

  if (nargin < 5)
    hard_thresh = 2e-2;
  end

  thresh = 1e-5;
  sub_coef = 2;
  ruffles = get_struct('ruffles', [1 nframes]); 

  for i=1:nframes

    pts = mymovie.(type).cortex(i).carth;

    if (isempty(pts))
      continue;
    end

    center = mymovie.(type).centers(:,i);
    axes_length = mymovie.(type).axes_length(:,i);
    orient = mymovie.(type).orientations(1,i);
  
    x = unique(pts(:,1));
    y = unique(pts(:,2));
    if (numel(x) < 5 | numel(y) < 5)
      continue;
    end

    conv = pts(convhull(pts(:,1),pts(:,2)),:);
    if (~ispolycw(conv(:,1),conv(:,2)))
      conv = conv([end:-1:1],:);
    end

    %figure;myplot(pts);hold on;
    %myplot(conv,'r');

    ell_pts = carth2elliptic(pts, center, axes_length, orient, 'radial');
    ell_conv = carth2elliptic(conv, center, axes_length, orient, 'radial');
    indx = find(ell_conv(1:end-1,1) > ell_conv(2:end,1));

    if (~isempty(indx))
      ell_conv = ell_conv([indx+1:end 1:indx],:);
    end

    tmp_conv = ell_conv([end 1:end 1],:);
    tmp_conv(1,1) = tmp_conv(1,1) - 2*pi;
    tmp_conv(end,1) = tmp_conv(end,1) + 2*pi;
    tmp_conv = interp1q(tmp_conv(:,1), tmp_conv(:,2), ell_pts(:,1));

    ell_conv = [ell_pts(:,1) tmp_conv];

    %figure;myplot(ell_pts);hold on;
    %myplot(ell_conv,'r');
      
    dist = ell_conv(:,2) - ell_pts(:,2);
    npts = length(dist);
    dist = [dist; dist(1:ceil(npts/2))];

    indxs = zeros(1,3);
    cands = zeros(1,2);
    local_max = 0;
    first = find(dist < thresh, 1);

    for j=first:length(dist)
      if (~isfinite(dist(j)))
        break;
      elseif (dist(j) < thresh & local_max > thresh)
        local_max = 0;
        indxs(end, 3) = j;

        if (indxs(1,2) == indxs(end,2)-npts)
          break;
        end

        indxs(end+1,:) = j;
      elseif (dist(j) > local_max)
        if (local_max < thresh)
          indxs(end, 1) = j-1;
        end

        local_max = dist(j);
        indxs(end, 2) = j;
      end
    end

    %keyboard
    
    indxs = indxs(1:end-1,:);

    max_dist = dist(indxs(:,2));
    %indxs = indxs(max_dist > mean(max_dist) + coef*std(max_dist),:);
    indxs = indxs(max_dist > hard_thresh,:);
    %bad_indxs = setdiff(indxs, good_indxs, 'rows');

    %ttest2(dist(good_indxs(:,2)),dist(bad_indxs(:,2)))

    %indxs = good_indxs;


    %%% SUBMAX

    %keyboard
    %figure;
    %plot([1:length(dist)], dist); hold on;

    new_indxs = zeros(0,3);
    
    for j=1:size(indxs,1)
      tmp = local_maxs(indxs(j,:), dist);
      tmp = tmp(1:2, tmp(3,:) > hard_thresh / sub_coef);

      nmaxs = size(tmp,2);
      if (nmaxs > 1)
        maxs = sort(tmp(1,:));
        mins = sort([tmp(2,:) indxs(j,3)]);

        for k=1:nmaxs
          new_indxs(end+1,:) = [mins(k) maxs(k) mins(k+1)];
        end
      else
        new_indxs(end+1,:) = indxs(j,:);
      end
      %scatter(tmp(1,:), dist(tmp(1,:)), 'k')
    end

    indxs = new_indxs - npts*(new_indxs > npts); 

    pts = [pts(indxs(:,2),:) pts(indxs(:,1),:) pts(indxs(:,3),:)];

    ruffles(i).carth = pts(:,1:2);
    ruffles(i).bounds = pts(:,3:end);
    ruffles(i).properties = dist(indxs(:,2));
  end

  mymovie.(type).ruffles = ruffles;

  %scatter(indxs(:,1), dist(indxs(:,1)), 'b')
  %scatter(indxs(:,2), dist(indxs(:,2)), 'r')
  %scatter(indxs(:,3), dist(indxs(:,3)), 'g')

  return;
end

function res = local_maxs(indxs, dist)

  first = indxs(1);
  final = indxs(3);

  dir = 1;
  prev = dist(first);
  mins = [];
  maxs = [];

  for i=first+1:final
    if (dist(i) > prev & dir == -1)
      mins = [mins i-1];
      dir = 1;
    elseif (dist(i) < prev & dir == 1)
      maxs = [maxs i-1];
      dir = -1;
    end
    prev = dist(i);
  end

  %scatter(maxs, dist(maxs), 'm')
  %scatter(mins, dist(mins), 'c')

  res = order_ruffles(dist, [first mins], maxs);  
  res = [res([2 1],:); diff(dist(res))];

  return;
end

function sort_mins = order_ruffles(path, mins, maxs)

  npts = size(mins,2);
  sort_mins = zeros(size([mins;maxs]));
  indexes = [];

  if (mins(1,1)<maxs(1,1))
    mins = mins(:,[2:end 1]);
  end

  mins_val = path(mins(1,:));
  maxs_val = path(maxs(1,:));

%  %keyboard
%
  d = -(bsxfun(@minus, mins_val, maxs_val'));
  d(isnan(d)) = min(d(:)) - 1;
  dindx = ones(npts);

  for i=1:npts

    found = false;
    mask = ones(npts);

    while (~found)
      td = d .* dindx .* mask;

      indxs = find(td==max(max(td)));

      if (td(indxs(1)) == 0)
        row = -1;
        break;
      end

      [rows,cols] = ind2sub([npts npts], indxs);

      jndx = 1;
      ncand = length(rows);
      if (ncand>1)
        dist = abs([rows - cols; rows - (cols+npts)]);
        jndx = find(dist==min(dist),1);

        if (jndx>ncand)
          jndx = jndx - ncand;
        end
      end

      row = rows(jndx);
      col = cols(jndx);

      found = valid_index(indexes, mins(1,row), maxs(1,col));

      if (~found)
        mask(row, col) = 0;
      end
    end

    if (row < 0)
      sort_mins = sort_mins(:, ~any(sort_mins == 0, 1));
      break;
    else
      sort_mins(:,i) = [mins(:,row);maxs(:,col)];
      indexes = [indexes [mins(1,row); maxs(1,col)]];
      dindx(row,:) = 0;
      dindx(:,col) = 0;
    end
  end
     

  return
%
%  % More efficient but buggy algorithm
%  uones = triu(dindx,1);
%  lones = tril(dindx);
%
%  minp = npts;
%  maxp = 1;
%  corner = 0;
%
%  %figure;
%  %plot(path,1:length(path),'g')
%  %hold on;
%
%  for i=1:npts
%
%    td = d .* dindx;
%    
%    if (~any(any(td>0)))
%      if (~any(any(dindx)))
%        break;
%      end
%      td = dindx;
%    end
%
%    indxs = find(td==max(max(td)));
%    [rows,cols] = ind2sub([npts npts], indxs);
%
%    jndx = 1;
%    ncand = length(rows);
%    if (ncand>1)
%      dist = abs([rows - cols; rows - (cols+npts)]);
%      jndx = find(dist==min(dist),1);
%
%      if (jndx>ncand)
%        jndx = jndx - ncand;
%      end
%    end
%
%    row = rows(jndx);
%    col = cols(jndx);
%
%    if (row < col)
%
%      if (row <= minp & col >= maxp)
%
%        dindx(1:minp-1,maxp+1:end) = 0;
%        dindx(1:end,1:row) = uones(1:end,1:row);
%        dindx(col:end,1:end) = uones(col:end,1:end);
%        dindx(col:end,1:row) = 1;
%
%        if (i==1)
%          dindx(row:col-1,row+1:col) = lones(row:col-1,row+1:col);
%        end
%
%        minp = row;
%        maxp = col;
%        corner = 1;
%      else
%        dindx(1:col-1,row+1:end) = lones(1:col-1,row+1:end);
%
%        if (row <= minp)
%          dindx(maxp:end,row+1:minp) = 0; 
%          minp = row;
%        elseif (col >= maxp)
%          dindx(maxp:col-1,1:minp) = 0; 
%          maxp = col;
%        end
%
%        if (corner~=1)
%          dindx(1:minp-1,maxp+1:end) = 1;
%        end
%      end
%    else
%
%      if (col <= minp & row >= maxp)
%      
%        dindx(maxp:end,1:minp) = 0;
%        dindx(1:end,row+1:end) = lones(1:end,row+1:end);
%        dindx(1:col-1,1:end) = lones(1:col-1,1:end);
%        dindx(1:col-1,row+1:end) = 1;
%
%        if (i==1)
%          dindx(col:row-1,col+1:row) = uones(col:row-1,col+1:row);
%        end
%
%        minp = col;
%        maxp = row;
%        corner = -1;
%      else
%        dindx(col:end,1:row) = uones(col:end,1:row);
%
%        if (col <= minp)
%          dindx(col:minp,maxp:end) = 0;
%          minp = col;
%        elseif (row >= maxp)
%          dindx(1:minp-1,maxp:row) = 0;
%          maxp = row;
%        end
%
%        if (corner~=-1)
%          dindx(maxp:end,1:minp) = 1;
%        end
%      end
%    end
%    %plot(mins_val(row),mins(1,row),'or')
%    %plot(maxs_val(col),maxs(1,col),'ob')
%
%    %keyboard
%
%    sort_mins(:,i) = [mins(:,row);d(row,col)];
%  end
%
%  return;
end

function ok = valid_index(pos, new_min, new_max)

  npts = size(pos,2);
  ok = false;

  if (npts == 0)
    ok = true;

    return;
  end

  min_indx = 2*npts + 1;
  max_indx = min_indx + 1;
  [values, indx] = sort([pos(1,:) pos(2,:) new_min new_max]);
  min_pos = find(indx == min_indx);

  if (min_pos == 1)
    ok = (indx(2) == max_indx && indx(3) <= npts) || ...
         (indx(end) == max_indx && indx(end-1) <= npts);
  elseif (min_pos == max_indx)
    ok = (indx(end-1) == max_indx && indx(end-2) <= npts) || ...
         (indx(1) == max_indx && indx(2) <= npts);
  else
    ok = (indx(min_pos-1) == max_indx && indx(min_pos+1) > npts) || ...
         (indx(min_pos+1) == max_indx && indx(min_pos-1) > npts);
  end

  return;
end
