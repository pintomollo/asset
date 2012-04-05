function [img, minval, maxval] = imnorm(img, minval, maxval, mode, minrange, maxrange)

  if (nargin == 1)
    minval = [];
    maxval = [];
    mode = '';
  elseif (nargin < 4)
    mode = '';
  end

  if (nargin < 6)
    minrange = 0;
    maxrange = 1;
  end

  if (isempty(img))
    return;
  end

  [h w] = size(img);

  if (isempty(minval))
    if (~isempty(mode) & mode(1) == 'r')
      minval = repmat(min(img,[],2), 1, w);
    elseif (~isempty(mode) & mode(1) == 'c')
      minval = repmat(min(img,[],1), h, 1);
    else
      minval = min(img(:));
    end
  end
  if (isempty(maxval))
    if (~isempty(mode) & mode(1) == 'r')
      maxval = repmat(max(img,[],2), 1, w);
    elseif (~isempty(mode) & mode(1) == 'c')
      maxval = repmat(max(img,[],1), h, 1);
    else
      maxval = max(img(:));
    end
  end

  img = (img - (minval - minrange)) .* ((maxrange - minrange) ./ (maxval - minval));
  img(img < minrange) = minrange;
  img(img > maxrange) = maxrange;

  return;
end
