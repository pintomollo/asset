function [ptsx, ptsy, indexes] = insert_ruffles(ptsx, ptsy, coords, path)

  if (nargin == 3)
    path = coords;
    coords = ptsy;

    ptsy = ptsx(:, 2);
    ptsx = ptsx(:, 1);
  end

  indexes = false(size(ptsx));

  for i = 1:size(coords, 1)
    if (~isempty(path{i}))
      insert = path{i}([1:end end-1:-1:1], :);

      indx = find(ptsx == coords(i, 1) & ptsy == coords(i, 2), 1);
      if (isempty(indx))
        dist = ((ptsx - coords(i,1)).^2 + (ptsy - coords(i, 2)).^2);
        [~, indx] = min(dist);
        indx = indx(1);
      end

      ptsx = [ptsx(1:indx); insert(:, 1); ptsx(indx:end)];
      ptsy = [ptsy(1:indx); insert(:, 2); ptsy(indx:end)];
      indexes = [indexes(1:indx-1); true(length(insert)+2, 1); indexes(indx+1:end)];
    end
  end

  if (nargout == 1)
    ptsx = [ptsx ptsy];
    ptsy = [];
  elseif (nargout == 2 & nargin == 3)
    ptsx = [ptsx ptsy];
    ptsy = indexes;
    indexes = [];
  end

  return;
end
