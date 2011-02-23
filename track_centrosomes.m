function mymovie = track_centrosomes(mymovie, opts)

  %hwait = waitbar(0,'Tracking centrosomes','Name','RECOS');

  if (opts.recompute | ~isfield(mymovie.data, 'eggshell') | isempty(mymovie.data.eggshell))
    
    mymovie = correct_dic_shift(mymovie, 'data', opts.segmentation_parameters.correction, opts);
  end

  [imgsize, nframes] = size_data(mymovie.data);

  sep_thresh = 5;
  weight = 0.45;
  
  centr1 = NaN(nframes, 2);
  centr2 = NaN(nframes, 2);

  tmp_pts = zeros(0, 3);

  for i=1:nframes

    nimg = i;
    
    img = load_data(mymovie.data, nimg);

    if (opts.verbosity == 3)
      imshow(img);
      drawnow
    end

    spots = detect_spots(img);

    switch (size(spots,1))
      case 0
        continue;
      case 1
        pts = spots(1,1:2);
        tmp_pts = [tmp_pts; pts i];
      otherwise
        pts = spots(1:2,1:2);
        tmp_pts = [tmp_pts; pts [i; i]];

        [centr1, centr2] = assign_spots(tmp_pts, centr1, centr2);
        tmp_pts = [centr1(i,: ) i; centr2(i,:) i];
    end

    if (opts.verbosity == 3)
      hold on;

      scatter(pts(:,2), pts(:,1))
    end
  end

  %if (isfield(mymovie, 'dic') & ~isempty(mymovie.dic))
  %  type = 'dic';
  %elseif (isfield(mymovie, 'markers') & ~isempty(mymovie.markers))
  %  type = 'markers';
  %end

  pts = [centr1(:,[2 1]) centr2(:, [2 1])];
  indx = find(~any(isnan(pts), 2), 1, 'last');

  try
  ell_pos = carth2elliptic(pts(indx,:), mymovie.data.centers(:,indx),mymovie.data.axes_length(:,indx),mymovie.data.orientations(1,indx));
  catch
    beep;beep;
    keyboard
  end
  if (abs(ell_pos(1)-pi) < abs(ell_pos(3)-pi))
    pts = pts(:,[3 4 1 2]);
  end

  for i=1:nframes
    mymovie.data.centrosomes(i).carth = [pts(i,1:2); pts(i,3:4)];
  end

  mymovie = resolve_dv(mymovie);
  % Done in resolve already
  %mymovie = carth2RECOS(mymovie);

  return;
end

function [centr1, centr2] = assign_spots(pts, centr1, centr2)

  npts = size(pts, 1);
  range = unique([min(pts(:,3)) max(pts(:,3))]);

  if (pts(1,3) == pts(2,3))
    
    mutual_dist = sqrt(bsxfun(@minus,pts(:,1),pts(:,1).').^2 + bsxfun(@minus,pts(:,2),pts(:,2).').^2);
    mutual_dist = mutual_dist + tril(Inf(size(mutual_dist)));
    mutual_dist = mutual_dist(1:end-2,3:end);

    indxs = [1 2];
    shift = 2;
  elseif (pts(end,3) == pts(end-1,3))
    pts = pts(end:-1:1,:);

    mutual_dist = sqrt(bsxfun(@minus,pts(:,1),pts(:,1).').^2 + bsxfun(@minus,pts(:,2),pts(:,2).').^2);
    mutual_dist = mutual_dist + tril(Inf(size(mutual_dist)));
    mutual_dist = mutual_dist(1:end-1,3:end);

    indxs = [2 1];
    shift = 2;
  else

    mutual_dist = sqrt(bsxfun(@minus,pts(:,1),pts(:,1).').^2 + bsxfun(@minus,pts(:,2),pts(:,2).').^2);
    mutual_dist = mutual_dist + tril(Inf(size(mutual_dist)));

    indxs = [1 0];
    shift = 0;
  end

  if (all(isnan(centr1(:)) & isnan(centr2(:))))
    real_indx = pts(indxs(1),3);
    centr1(real_indx,:) = pts(indxs(1),1:2);
    if (indxs(2) ~= 0)
      real_indx = pts(indxs(1),3);
      centr2(real_indx,:) = pts(indxs(2),1:2);
    end
  end

  assign = munkres(mutual_dist);
  for i=1:length(assign)
    if (assign(i) == 0)
      continue;
    end

    new_indx = assign(i) + shift;
    real_indx = pts(new_indx,3);
    if (i == indxs(1))
      centr1(real_indx,:) = pts(new_indx,1:2);
      indxs(1) = new_indx;
    else
      centr2(real_indx,:) = pts(new_indx,1:2);
      indxs(2) = new_indx;
    end
  end

  %keyboard

  centr1 = fill_gaps(centr1, range); 
  centr2 = fill_gaps(centr2, range); 

  return;
end

function [pos] = fill_gaps(pos, indxs)

  if (length(indxs) == 1)
    return;
  end

  indxs = [indxs(1):indxs(2)];
  indxs = indxs(~isnan(pos(indxs,1)));
  
  npts = length(indxs);

  dist = diff(indxs);  
  for i=1:npts-1
    if (dist(i) > 1)
      ninter = dist(i)-1;
      dmov = pos(indxs(i+1),:) - pos(indxs(i),:);
      dx = [1:ninter].';
      new_pts = pos(indxs(i)*ones(1,ninter),:) + dx * dmov / dist(i);
      pos(indxs(i)+1:indxs(i+1)-1,:) = new_pts; 
    end
  end

  return;
end
