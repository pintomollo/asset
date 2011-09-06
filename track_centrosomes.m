function mymovie = track_centrosomes(mymovie, opts)
% TRACK_CENTROSOMES detects and tracks over time the position of the two centrosomes
% in a fluorescence time-lapse recording. This algorithm is designed to find two and
% only two "particles".
%
%   MYMOVIE = TRACK_CENTROSOMES(MYMOVIE, OPTS) detects particles in the "data"
%   channel of MYMOVIE, identifies the centrosomes and tracks them over time. These
%   results are stored in the "centrosomes" field of "data".
%
%   [...] = TRACK_CENTROSOMES(MYMOVIE) uses the default value of OPTS as provided by
%   get_struct('ASSET').
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 01.09.2011

  % Check whether we actually do have to compute anything 
  if (isfield(mymovie.data, 'centrosomes') & ~isempty(mymovie.data.centrosomes) & ~opts.recompute)
    return;
  end

  % The the recording size
  [nframes, imgsize] = size_data(mymovie.data);

  % Identify the candidate spots
  spots = detect_spots(mymovie.data, opts);

  % Start a very simple tracking, base don the fact that we look for two and only
  % two particles

  % Initialize some variables which should not be hard-coded here !!
  sep_thresh = 5;
  weight = 0.45;

  % Assign variables for a slight speed-up
  centr1 = NaN(nframes, 2);
  centr2 = NaN(nframes, 2);

  % We'll start by collecting the position of the centrosomes frame by frame
  tmp_pts = zeros(0, 3);
  for i=1:nframes

    % Extract the current set of candidates
    nimg = i;
    tmp_spots = spots{i};

    % Group them happily ! 
    switch (size(tmp_spots,1))
      % Nothing here...
      case 0
        continue;

      % We cannot decide which centrosome we see until we have located both of them
      % in a single frame
      case 1
        pts = tmp_spots(1,1:2);
        tmp_pts = [tmp_pts; pts i];

      % We have both in the current frame !
      otherwise
        % Keep only the two best candidates (ordered by detect_spots.m already)
        pts = tmp_spots(1:2,1:2);
        tmp_pts = [tmp_pts; pts [i; i]];

        % Re-organize all the previous particles into two paths
        [centr1, centr2] = assign_spots(tmp_pts, centr1, centr2);
        % Store the latest time point for further tracking
        tmp_pts = [centr1(i,: ) i; centr2(i,:) i];
    end

    % Amazingly fancy display
    if (opts.verbosity == 3)
      imshow(imnorm(double(load_data(mymovie.data, nimg))));
      hold on;
      scatter(pts(:,2), pts(:,1))
      hold off;
    end
  end

  % Group the two centrosomes and resolve the X-Y vs I-J problem of indexes
  pts = [centr1(:,[2 1]) centr2(:, [2 1])];
  % Keep the last frame in which we have both
  indx = find(~any(isnan(pts), 2), 1, 'last');

  % Get the elliptical position of both centrosomes to identify which one is the
  % anterior one
  ell_pos = carth2elliptic([pts(indx,1:2); pts(indx, 3:4)], mymovie.data.centers(:,indx),mymovie.data.axes_length(:,indx),mymovie.data.orientations(1,indx));
  % The anterior one is closer to pi but should be in 2nd position
  if (abs(ell_pos(1, 1)-pi) < abs(ell_pos(2, 1)-pi))
    pts = pts(:,[3 4 1 2]);
  end

  % Store the results in mymovie
  for i=1:nframes
    mymovie.data.centrosomes(i).carth = [pts(i,1:2); pts(i,3:4)];
  end

  % Check whether there is need to invert the D-V axis
  mymovie = resolve_dv(mymovie);

  return;
end

% This function tracks two centrosomes over several frames, briging possible gaps.
% The new positions are in pts, the previously assigned ones in centrN
function [centr1, centr2] = assign_spots(pts, centr1, centr2)

  % The number of points to assign
  npts = size(pts, 1);
  % The range of frames in which we play
  range = unique([min(pts(:,3)) max(pts(:,3))]);

  % If we know the position of both centrosomes in the first frame, that's easy !
  if (pts(1,3) == pts(2,3))
    
    % Simply compute the all-to-all distance
    mutual_dist = sqrt(bsxfun(@minus,pts(:,1),pts(:,1).').^2 + bsxfun(@minus,pts(:,2),pts(:,2).').^2);
    % Restrict the possible assignations
    mutual_dist = mutual_dist + tril(Inf(size(mutual_dist)));
    mutual_dist = mutual_dist(1:end-2,3:end);

    % Two variables necessary for the tracking. indxs indicates which row in
    % mutual_dist belongs to which track while shift indicates how many points we
    % removed from the assignation.
    indxs = [1 2];
    shift = 2;

  % If we know both at the end, we can back-track them
  elseif (pts(end,3) == pts(end-1,3))
    pts = pts(end:-1:1,:);

    mutual_dist = sqrt(bsxfun(@minus,pts(:,1),pts(:,1).').^2 + bsxfun(@minus,pts(:,2),pts(:,2).').^2);
    mutual_dist = mutual_dist + tril(Inf(size(mutual_dist)));
    mutual_dist = mutual_dist(1:end-1,3:end);

    indxs = [2 1];
    shift = 2;

  % Otherwise, we simply assign as we can
  else
    mutual_dist = sqrt(bsxfun(@minus,pts(:,1),pts(:,1).').^2 + bsxfun(@minus,pts(:,2),pts(:,2).').^2);
    mutual_dist = mutual_dist + tril(Inf(size(mutual_dist)));

    indxs = [1 0];
    shift = 0;
  end

  % In this case, we have no information yet about the centrosomes.
  % So we need to choose arbitrarly the anterior and posterior ones.
  if (all(isnan(centr1(:)) & isnan(centr2(:))))
    % Get the actual index of the first centrosome
    real_indx = pts(indxs(1),3);
    % Store it as a reference
    centr1(real_indx,:) = pts(indxs(1),1:2);

    % If we have a second candidate, store it as well
    if (indxs(2) ~= 0)
      real_indx = pts(indxs(1),3);
      centr2(real_indx,:) = pts(indxs(2),1:2);
    end
  end

  % Perform the minimum distance assignation
  assign = munkres(mutual_dist);

  % Now track the centrosomes !
  for i=1:length(assign)
    if (assign(i) == 0)
      continue;
    end

    % Get the actual index, after the shift due to the inital positions that we removed
    new_indx = assign(i) + shift;
    real_indx = pts(new_indx,3);

    % Assign them to the correct track
    if (i == indxs(1))
      centr1(real_indx,:) = pts(new_indx,1:2);

      % Store which is the next index in this track
      indxs(1) = new_indx;
    else
      centr2(real_indx,:) = pts(new_indx,1:2);
      indxs(2) = new_indx;
    end
  end

  % Close possible gaps using simple linear interpolation
  centr1 = fill_gaps(centr1, range); 
  centr2 = fill_gaps(centr2, range); 

  return;
end

% Here we bridge the position of centrosomes over frames in which they are missing
function [pos] = fill_gaps(pos, indxs)

  if (length(indxs) == 1)
    return;
  end

  % Get the indexes in which we do have the positions
  indxs = [indxs(1):indxs(2)];
  indxs = indxs(~isnan(pos(indxs,1)));
  
  % Get the number of frames
  npts = length(indxs);

  % Get the dT between each frame
  dist = diff(indxs);  
  for i=1:npts-1
    % If there is indeed a gap, bridge it using a simple linear interpolation
    if (dist(i) > 1)
      ninter = dist(i)-1;
      dmov = pos(indxs(i+1),:) - pos(indxs(i),:);
      dx = [1:ninter].';
      new_pts = pos(indxs(i)*ones(1,ninter),:) + dx * dmov / dist(i);

      % Store the results
      pos(indxs(i)+1:indxs(i+1)-1,:) = new_pts; 
    end
  end

  return;
end
