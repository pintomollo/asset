function [values, perp_path, dpos] = perpendicular_sampling(img, path, perp_path, dpos, opts)

  if (nargin == 2)
    [path, opts] = deal(img, path);
    img = [];
    perp_path = [];
    dpos = [];
  elseif (nargin == 3)
    opts = perp_path;
    perp_path = [];
    dpos = [];
  elseif (nargin == 4)
    if (all(size(path) == size(perp_path)))
      opts = dpos;
      dpos = [];
    else
      opts = dpos;
      dpos = perp_path;
      perp_path = [];
    end
  end

  if (isempty(perp_path))

    [dist, tot_dist] = carth2linear(path);
    dist = dist * tot_dist;

    deriv = differentiator([dist(:), dist(:)], path, 'circular');
    perp_path = [-deriv(:,2), deriv(:,1)];

    nulls = all(perp_path == 0, 2);

    if (any(nulls))
      dpts = diff(path([1:end 1], :));
      perp_path(nulls, :) = dpts(nulls, :);
    end

    perp_path = bsxfun(@rdivide, perp_path, hypot(perp_path(:,1), perp_path(:,2))); 
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

  if (isempty(img))
    values = perp_path;
    perp_path = dpos;

    dpos = [];

    return;
  end

  all_pos_x = bsxfun(@plus, perp_path(:,1) * dpos, path(:,1));
  all_pos_y = bsxfun(@plus, perp_path(:,2) * dpos, path(:,2));

  all_pos_x = all_pos_x(:);
  all_pos_y = all_pos_y(:);

  values = bilinear(img, all_pos_x, all_pos_y);
  values = reshape(values, size(path,1), []);

  return;
end
