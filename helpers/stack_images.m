function [window, params] = stack_images(varargin)

  imgs = cell(0,1);
  for i=1:nargin-1
    if (iscell(varargin{i}))
      imgs = [imgs; varargin{i}];
    else
      imgs{end+1,1} = varargin{i};
    end
  end

  params = varargin{end};
  if (~isnumeric(params) || numel(params)~=length(imgs))
    if (iscell(params))
      imgs = [imgs; params];
    else
      imgs{end+1,1} = params;
    end

    params = ones(length(imgs), 1);
  end

  nimgs = length(imgs);

  sizes = NaN(nimgs, 2);
  for i=1:nimgs
    sizes(i,:) = size(imgs{i});
  end

  after = sizes(:,1) - params;
  before = max(params);
  shift = before - (params - 1);

  window_size = [before+max(after) max(sizes(:,2))];
  window = NaN([window_size nimgs]);

  dw = round((window_size(2) - sizes(:,2))/2);
  for i=1:nimgs
    window(shift(i):before+after(i,1),dw(i)+1:dw(i)+sizes(i,2),i) = imgs{i};
  end

  params = params - min(params) + 1;

  return;
end
