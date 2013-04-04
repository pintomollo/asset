function [window, shift, pos] = stack_images(varargin)

  params = [];
  thresh = 0;
  resolution = 1;
  imgs = cell(0,1);
  for i=1:nargin
    if (iscell(varargin{i}))
      imgs = [imgs; varargin{i}];
    elseif (numel(varargin{i})==1)
      if (thresh == 0)
        thresh = varargin{i};
      else
        resolution = varargin{i};
      end
    else
      [n,m,o] = size(varargin{i});

      if ((n==1 || m==1) && o==1)
        params = varargin{i};
      else
        imgs{end+1,1} = varargin{i};
      end
    end
  end

  if (isempty(params) || any(~isnumeric(params)) || numel(params)~=length(imgs))
    params = ones(length(imgs), 1);
  end
  params = params(:);

  nimgs = length(imgs);

  sizes = NaN(nimgs, 3);
  for i=1:nimgs
    sizes(i,1) = size(imgs{i}, 1);
    sizes(i,2) = size(imgs{i}, 2);
    sizes(i,3) = size(imgs{i}, 3);
  end

  after = sizes(:,1) - params;
  before = max(params);
  shift = before - (params - 1);

  window_size = [before+max(after) max(sizes(:,2))];
  window = NaN([window_size sum(sizes(:,3))]);

  dw = round((window_size(2) - sizes(:,2))/2);
  count = 1;
  for i=1:nimgs
    window(shift(i):before+after(i,1),dw(i)+1:dw(i)+sizes(i,2),count:count+sizes(i,3)-1) = imgs{i};
    count = count + sizes(i,3);
  end

  params = params - min(params) + 1;

  if (thresh > 0)
    half = ((size(window, 2)-1) / 2);

    pos = [-half:half]*resolution;

    counts = nanmean(isfinite(window), 3);
    bads = (counts <= thresh);

    rows = ~all(bads, 2);
    first = find(rows, 1, 'first');
    last = find(rows, 1, 'last');
    window = window(first:last, :, :);

    bads = bads(first:last, :);
    cols = ~all(bads, 1);
    left = find(cols, 1, 'first');
    right = size(bads,2) - find(cols, 1, 'last') + 1;
    width = max(left, right);

    window = window(:, width:end-width+1, :);
    pos = pos(width:end-width+1);
  end

  return;
end
