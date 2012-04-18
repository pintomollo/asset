function [ptsx, ptsy, indexes] = insert_ruffles(ptsx, ptsy, path)

  if (nargin == 2)
    path = ptsy;

    ptsy = ptsx(:, 2);
    ptsx = ptsx(:, 1);
  end

  indexes = false(size(ptsx));

  for i = 1:length(path)
    if (~isempty(path{i}))
      insert = path{i}([1:end end-1:-1:1], :);

      indx = find(ptsx == insert(1, 1) & ptsy == insert(1, 2), 1);
      if (isempty(indx))
        dist = ((ptsx - insert(1, 1)).^2 + (ptsy - insert(1, 2)).^2);
        [junk, indx] = min(dist);
        indx = indx(1);

        if (dist > 1e-10)
          insert_pt = [ptsx(indx) ptsy(indx)];
          insert = [insert_pt; insert; insert_pt];
        end
      end

      ptsx = [ptsx(1:indx-1); insert(:, 1); ptsx(indx+1:end)];
      ptsy = [ptsy(1:indx-1); insert(:, 2); ptsy(indx+1:end)];
      indexes = [indexes(1:indx-1); true(size(insert, 1), 1); indexes(indx+1:end)];
    end
  end

  if (nargout == 1)
    ptsx = [ptsx ptsy];
    ptsy = [];
  elseif (nargout == 2 && nargin == 2)
    ptsx = [ptsx ptsy];
    ptsy = indexes;
    indexes = [];
  end

  return;
end
