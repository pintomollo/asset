function [domain, pos, center_indx] = align_domain(mymovie, opts, path)

  if (nargin == 3)
    tmp = path;
    path = opts;
    opts = tmp;
  else
    path = mymovie.data.domain;
  end

  [img, pos] = gather_quantification(mymovie, opts);
  %path = path * opts.pixel_size;

  domain = [];
  sides = zeros(1, 2);
  [h,w] = size(img);
  resolution = median(diff(pos));
  rel_pos = pos;

  center_indx = find(pos == 0);
  %center = (path(:, 1) - center_indx);
  center = -1;

  %keyboard

  for i=h:-1:1
    line = img(i,:);
    if (~isnan(path(i,1)))
      center = round(path(i, 1));
    elseif (center < 0)
      center = center_indx;
    end

    min_indx = find(~isnan(line), 1, 'first');
    max_indx = find(~isnan(line), 1, 'last');

    min_pos = (min_indx - center) * resolution;
    max_pos = (max_indx - center) * resolution;

    if (min_pos < sides(1))
      nmin = round((sides(1) - min_pos) / resolution);
      sides(1) = min_pos;
    else
      nmin = 0;
    end
    if (max_pos > sides(2))
      nmax = round((max_pos - sides(2)) / resolution);
      sides(2) = max_pos;
    else
      nmax = 0;
    end     
    if (i==h)
      domain = img(i, min_indx:max_indx);
    else
      nrows = size(domain, 1);
      domain = [NaN(nrows, nmin), domain, NaN(nrows, nmax)];

      first_indx = round((min_pos - sides(1))/resolution) + 1;
      %last_indx = (sides(2) - max_pos) / resolution;

      domain(end+1, first_indx:first_indx+(max_indx-min_indx)) = line(min_indx:max_indx);

    end
  end
  domain = flipud(domain);
  pos = [fliplr([0:-resolution:sides(1)]) resolution:resolution:(sides(2)+1)];
  pos = pos(1:size(domain, 2));

  return;
end
