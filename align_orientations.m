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

  if (nargin < 3)
    if (ndims(orients) <= 2 && size(orients,2) > 1)
      orients = orients.';
    end

    dim = 1;
  end
  if (nargin < 2 | isempty(aim))
    aim = NaN;
  end
  last_dim = ndims(orients) + 1;
  if (last_dim == 3 && size(orients, 2) == 1)
      last_dim = 2;
  end

  % Test
  
  [m, n, o, p, q] = size(orients);
  orient_size = [m, n, o, p, q];

  rep_size = ones(1, last_dim);
  rep_size(dim) = orient_size(dim);
  rep_size(last_dim) = 3;

  if (isnan(aim))
    avg = median(orients, dim);
  else
    avg = aim;
  end
  orients = cat(last_dim, orients, (orients - pi), (orients + pi));
  dist = abs(orients - repmat(avg, rep_size));

  orient_size(last_dim) = 3;

  [junk, indx] = min(dist,[],last_dim);
  str = 'indx = sub2ind(size(orients)';
  for i=1:last_dim-1
        tmp_size = orient_size;
        tmp_size(i) = 1;
        tmp_size(last_dim) = 1;
      
        vect_size = ones(1,last_dim);
        vect_size(i) = orient_size(i);
      
        dim_indx = reshape(1:orient_size(i), vect_size);
        dim_indx = repmat(dim_indx, tmp_size);
        
      str = [str ', [' num2str(dim_indx(:).') ']'];
  end
  
  eval([str ', [' num2str(indx(:).') ']);']);
  
  orient_size(last_dim) = 1;
  orients = reshape(orients(indx), orig_size);

  if (~isempty(mystruct))
    mystruct.orientations = orients;
    orients = mystruct;
  end

  return;
end
