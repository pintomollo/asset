function h = imagesp(img, percent, varargin)

  if (nargin == 1)
    percent = 99;
    bounds = prctile(img(:), [100 - percent, percent]);
  elseif (nargin > 2)
    h = imagesc(img, percent, varargin{:});

    return;
  elseif (numel(percent) == 2)
    bounds = percent;
  else
    bounds = prctile(img(:), [100 - percent, percent]);
  end

  hh = imagesc(img, bounds);
  if (nargout == 1)
    h = hh;
  end

  return;
end
