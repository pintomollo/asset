function orients = align_orientations(orients, aim, dim)

  orig_size = size(orients);

  if (isstruct(orients))
    if (isfield(orients, 'orientations'))
      mystruct = orients;
      orients = orients.orientations;
      orig_size = size(orients);
    else
      return;
    end
  else
    mystruct = [];
  end

  if (nargin == 1)
    dim = 1;
    aim = [];
  elseif (nargin == 2)
    dim = 1;
  end

  orient_size = size(orients);
  if (dim == 1 & orient_size(dim) == 1)
    dim = find(orient_size > 1, 1);

    if (isempty(dim))
      return;
    end
  elseif (dim > numel(orient_size) | orient_size(dim) == 1)
    return;
  end

  perm_dim = [1:length(orient_size)];
  perm_dim(dim) = 1;
  perm_dim(1) = dim;

  orients = permute(orients, perm_dim);
  orients = reshape(orients, orient_size(dim), []);
  valids = isfinite(orients);
  orients(~valids) = 0;

  if (isfinite(aim))
    if (numel(aim) ~= 1)
      aim = permute(aim, perm_dim);
      aim = reshape(aim, orient_size(dim), []);
      aim = [aim;, repmat(aim(1,:), orient_size(dim) - size(aim, 1), 1)];
    end
    aim = aim(valids, :);
  else
    sorts = sort(orients, 1);
    aim = sorts(ceil(orient_size(dim)/2), :);
  end

  orients = cat(3, orients, (orients - pi), (orients + pi));
  dist = abs(bsxfun(@minus, orients, aim));

  sizes = size(dist);

  [junk, mins] = min(dist, [], 3);
  i = repmat([1:sizes(1)].', 1, sizes(2));
  j = repmat([1:sizes(2)], sizes(1), 1);
  indexes = sub2ind(sizes, i, j, mins);

  orients = orients(indexes);
  orients(~valids) = NaN;

  orients = reshape(orients, orient_size(perm_dim));
  orients = ipermute(orients, perm_dim);

  if (~isempty(mystruct))
    mystruct.orientations = orients;
    orients = mystruct;
  end

  return;
end
