function [domain, pos, center_indx] = align_domain(mymovie, opts, path)

  if (nargin == 3)
    tmp = path;
    path = opts;
    opts = tmp;
  else
    path = mymovie.data.domain;
  end

  if (isstruct(mymovie))
    [img, pos] = gather_quantification(mymovie, opts);
  else
    img = mymovie;
    pos = [1:size(img, 2)];
  end
  %path = path * opts.pixel_size;

  [h,w] = size(img);

  valids = ~isnan(img);
  path = round(path(:,1));
  last = find(~isnan(path), 1, 'last');
  path(last+1:end) = path(last);
  first = find(~isnan(path), 1, 'first');
  path(1:first-1) = path(first);
  width = NaN(h, 2);;

  indxs = [1:w];
  for i=1:h
    tmp = indxs(1:path(i));
    left = min(tmp(valids(i, 1:path(i))));
    tmp = indxs(path(i):end);
    right = max(tmp(valids(i, path(i):end)));

    if (isempty(left))
      left = path(i);
    end
    if (isempty(right))
      right = path(i);
    end

    width(i,:) = [left, right];
  end
  shift = bsxfun(@minus, width, path);
  center_indx = max(abs(shift(:))) + 1;

  domain = NaN(h, 2*center_indx-1);
  resolution = median(diff(pos));
  pos = [-(center_indx-1):center_indx-1] * resolution;
  %rel_pos = pos;

  for i=1:h
    indxs = [width(i,1) : width(i,2)];
    domain(i,indxs - path(i) + center_indx) = img(i, indxs);
  end

  %center_indx = find(pos == 0);
  %center = (path(:, 1) - center_indx);
  %center = -1;

  %keyboard

  %for i=h:-1:1
  %  line = img(i,:);
  %  if (~isnan(path(i,1)))
  %    center = round(path(i, 1));
  %  elseif (center < 0)
  %    center = center_indx;
  %  end

  %  min_indx = find(~isnan(line), 1, 'first');
  %  max_indx = find(~isnan(line), 1, 'last');

  %  min_pos = (min_indx - center) * resolution;
  %  max_pos = (max_indx - center) * resolution;
%
  %  if (min_pos < sides(1))
  %    nmin = round((sides(1) - min_pos) / resolution);
  %    sides(1) = min_pos;
  %  else
  %    nmin = 0;
  %  end
  %  if (max_pos > sides(2))
  %    nmax = round((max_pos - sides(2)) / resolution);
  %    sides(2) = max_pos;
  %  else
  %    nmax = 0;
  %  end     
  %  if (i==h)
  %    domain = img(i, min_indx:max_indx);
  %  else
  %    nrows = size(domain, 1);
  %    domain = [NaN(nrows, nmin), domain, NaN(nrows, nmax)];

  %    first_indx = round((min_pos - sides(1))/resolution) + 1;
      %last_indx = (sides(2) - max_pos) / resolution;
%
  %    domain(end+1, first_indx:first_indx+(max_indx-min_indx)) = line(min_indx:max_indx);

  %  end
  %end
  %domain = flipud(domain);
  %pos = [fliplr([0:-resolution:sides(1)]) resolution:resolution:(sides(2)+1)];
  %pos = pos(1:size(domain, 2));

  return;
end
