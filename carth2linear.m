function [linear, total_distance] = carth2linear(pts_x, pts_y, ruffles_path)

  if (nargin == 1)
    pts_y = [];
  end
  if (nargin < 3)
    ruffles_path = {};
  end

  if (iscell(pts_y))
    ruffles_path = pts_y;
    pts_y = [];
  end

  if (isempty(pts_y))
    pts_y = pts_x(:, 2);
    pts_x = pts_x(:, 1);
  end

  circular = true;

  [pts_x, pts_y] = insert_ruffles(pts_x, pts_y, ruffles_path);

  if (pts_x(1) ~= pts_x(end) | pts_y(1) ~= pts_y(end))
    pts_x = [pts_x; pts_x(1)];
    pts_y = [pts_y; pts_y(1)];

    circular = false;
  end

  distance = sqrt(diff(pts_x).^2 + diff(pts_y).^2);
  cum_dist = [0; cumsum(distance)];

  total_distance = cum_dist(end);

  if (total_distance == 0)
    linear = 0;

    return;
  end

  linear = cum_dist / total_distance;

  if (~circular)
    linear = linear(1:end-1);
  end

  return;
end
