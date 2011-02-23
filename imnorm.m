function [img, minval, maxval] = imnorm(img, minval, maxval, mode)

  if (nargin == 1)
    minval = [];
    maxval = [];
    mode = '';
  elseif (nargin < 4)
    mode = '';
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
  
  img = (img - minval) ./ (maxval - minval);
  img(img<0) = 0;
  img(img>1) = 1;

  return;
end
