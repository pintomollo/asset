function path = interp_elliptic(pts, pos)

  %keyboard

  dangle = mean(diff(pos(:,1))); 
  min_angle = min(pos(:,1)) - dangle;
  max_angle = max(pos(:,1)) + dangle;

  pts(pts(:,1) < min_angle,1) = pts(pts(:,1) < min_angle,1) + 2*pi;
  pts(pts(:,1) > max_angle,1) = pts(pts(:,1) > max_angle,1) - 2*pi;

  indx = find(pts(1:end-1,1) > pts(2:end,1), 1);
  if (~isempty(indx))
    pts = pts([indx+1:end 1:indx],:);
  end

  [npts, ndim] = size(pts);
  start = find(pts(:,1) - 2*pi < pos(1,1), 1, 'last');
  ends = find(pts(:,1) + 2*pi > pos(end,1), 1, 'first');

  pts = pts([start:end 1:end 1:ends],:);
  pts(1:npts-start+1,1) = pts(1:npts-start+1,1) - 2*pi;
  pts(end-ends+1:end,1) = pts(end-ends+1:end,1) + 2*pi;

  path = interp1q(pts(:,1), pts(:,2:end), pos);
  path = [pos path];

  return;
end
