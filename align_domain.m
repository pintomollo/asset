function [domain, ruffles, pos, center_indx] = align_domain(mymovie, opts, path)

  if (nargin == 3)
    tmp = path;
    path = opts;
    opts = tmp;
  else
    path = mymovie.data.domain;
  end

  if (isstruct(mymovie))
    [img, ruffles, pos] = gather_quantification(mymovie, opts);
  else
    img = mymovie;
    pos = [1:size(img, 2)];
  end
  %path = path * opts.pixel_size;

  [h,w] = size(img);

  valids = ~isnan(img);
  path = round(path(:,1));

  if (any(isnan(path(:))))
    
    slope = diff(path);
    x0 = find(~isnan(slope), 1, 'first');
    s0 = slope(x0);

    if (s0 == 0)
      path(1:x0-1) = path(x0);
    else
      plateau = find(slope == 0, 1, 'first');

      fact = 0.5;
      if (isempty(plateau))
        plateau = h;
      end

      C = path(plateau);
      A = fact*2*(path(x0) - C);
      B = s0 * 4 / A;

      sig = A./(1 + exp(-B*([1:x0-1] - x0))) + C + A*((1-fact)/(fact*2));
      path(1:x0-1) = round(sig);
    end

    %figure;plot(path, 1:h, 'b');
    %figure;imagesc(img);
    %hold on;plot(path, 1:h, 'r')

    last = find(~isnan(path), 1, 'last');
    path(last+1:end) = path(last);
    %first = find(~isnan(path), 1, 'first');
    %path(1:first-1) = path(first);
  end

  center_indx = ceil(w/2);
  domain = NaN(h, w);

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

    npts = right - left + 1;
    nleft = min(round(npts/2) + 1, center_indx);
    nright = npts - nleft;

    nfirst = max(path(i) - nleft + 1, left);
    nlast = min(path(i) + nright, right);

    domain(i, center_indx - nleft + 1 : center_indx + nright) = [img(i, nlast+1:right) img(i, nfirst:nlast), img(i, left:nfirst-1)];
  end

  resolution = median(diff(pos));
  pos = ([1:w] - center_indx) * resolution;

  return;
end
