function [perps_vals, perps_paths, dpos] = perpendicular_sampling(img, path, dpos)

  if (nargin < 3)
    dpos = [];
  end

  if (isempty(dpos))
    bin_dist = 64;
    bin_step = 1;
    dpos = [-bin_dist:bin_step:bin_dist];
  elseif (numel(dpos) == 2 & dpos(1) > dpos(2))
    bin_dist = dpos(1);
    bin_step = dpos(2);
    dpos = [-bin_dist:bin_step:bin_dist];
  end

  gaps = all(isnan(path), 2);
  dist = sqrt(diff(path(:,1)).^2 + diff(path(:,2)).^2);

  dist(isnan(dist)) = 0;

  cum_dist = [0; cumsum(dist)];

  stretchs = diff([0; cum_dist(gaps)]);

  dist(gaps(1:end-1)) = -stretchs(1:end-1);
  cum_dist = [0; cumsum(dist)];

  stretchs = cum_dist(gaps);

  cum_dist(gaps) = NaN;
  cum_dist([true; gaps(1:end-1)]) = 0;

  stretchs = find(gaps);
  nstretchs = length(stretchs);

  if (nstretchs == 0)
    nstretchs = 1;
    stretchs = size(path, 1) + 1;
  end

  perps_vals = cell(nstretchs, 1);
  perps_paths = cell(nstretchs, 1);
  prev = 1;

  for i=1:nstretchs
    indxs = [prev:stretchs(i)-1];
    curr_path = path(indxs,:);
    curr_dist = cum_dist(indxs);

    [curr_dist, indxs] = unique(curr_dist);
    curr_path = curr_path(indxs,:);

    dist = [0:floor(curr_dist(end))];
    curr_path = interp1(curr_dist.', curr_path, dist);

    dist = sqrt(diff(curr_path(:,1)).^2 + diff(curr_path(:,2)).^2);
    curr_dist = [0; cumsum(dist)];

    dist = [0:floor(curr_dist(end))];
    curr_path = interp1(curr_dist.', curr_path, dist);

    deriv = differentiator([dist(:), dist(:)], curr_path);
    perp_path = [-deriv(:,2), deriv(:,1)];

    nulls = all(perp_path == 0, 2);

    if (any(nulls))
      dpts = diff(curr_path([1:end 1], :));
      perp_path(nulls, :) = dpts(nulls, :);
    end

    perp_path = bsxfun(@rdivide, perp_path, hypot(perp_path(:,1), perp_path(:,2))); 

    all_pos_x = bsxfun(@plus, perp_path(:,1) * dpos, curr_path(:,1));
    all_pos_y = bsxfun(@plus, perp_path(:,2) * dpos, curr_path(:,2));

    all_pos_x = all_pos_x(:);
    all_pos_y = all_pos_y(:);

    values = bilinear_mex(double(img), all_pos_x, all_pos_y);
    values = reshape(values, size(curr_path,1), []);

    perps_vals{i} = values;
    perps_paths{i} = [all_pos_x, all_pos_y];

    prev = stretchs(i)+1;
  end

  if (length(perps_vals) == 1)
    perps_vals = perps_vals{1};
    perps_paths = perps_paths{1};
  end

  return;
end
