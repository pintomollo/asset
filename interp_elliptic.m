function [pos, path] = interp_elliptic(pts_t, pts_r, pos)
% LOAD_TRACKINGS loads the content of the .shapes files into the provided structure.
%
%   TRACKINGS = LOAD_TRACKINGS(TRACKINGS, OPTS) returns the same TRACKING structure
%   to which the content of the .shapes files will be added.
%
%   [TRACKINGS, NFRAMES] = LOAD_TRACKINGS(...) also returns the maximum number of
%   frames among all the .shapes files.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 23.03.2011

  %keyboard
  if (nargin == 2)
    pos = pts_r;
    pts_r = pts_t(:, 2:end);
    pts_t = pts_t(:, 1);
  end

  indx = find(pts_t(1:end-1) > pts_t(2:end), 1);
  if (~isempty(indx))
    pts_t = pts_t([indx+1:end 1:indx]);
    pts_r = pts_r([indx+1:end 1:indx], :);
  end

  if (numel(pos) == 1)
    resample = pos;
    pos = pts_t;
    for i=1:resample
      tmp_pos = NaN(size(pos,1)*2, 1);
      first = pos(1);
      first = first + 2*pi;

      mean_pos = mean([pos,[pos(2:end,:); first]], 2);
      tmp_pos(1:2:end-1,:) = pos;
      tmp_pos(2:2:end,:) = mean_pos;

      pos = tmp_pos;
    end
  end

  dangle = mean(diff(pos(:,1))); 
  min_angle = min(pos(:,1)) - dangle;
  max_angle = max(pos(:,1)) + dangle;

  pts_t(pts_t < min_angle) = pts_t(pts_t < min_angle) + 2*pi;
  pts_t(pts_t > max_angle) = pts_t(pts_t > max_angle) - 2*pi;

  npts = length(pts_t);
  start = find(pts_t - 2*pi < pos(1,1), 1, 'last');
  ends = find(pts_t + 2*pi > pos(end,1), 1, 'first');

  pts_t = pts_t([start:end 1:end 1:ends]);
  pts_r = pts_r([start:end 1:end 1:ends], :);
  pts_t(1:npts-start+1) = pts_t(1:npts-start+1) - 2*pi;
  pts_t(end-ends+1:end) = pts_t(end-ends+1:end) + 2*pi;

  path = interp1q(pts_t, pts_r, pos);

  if (nargout == 1)
    pos = [pos path];
    path = [];
  end

  return;
end
