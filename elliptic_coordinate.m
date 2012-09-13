function [elliptic_img] = elliptic_coordinate(img, varargin)

  [center, axes_length, orient, safety, ellsize, circular, type] = parse_inputs(varargin{:});
  if (any(~isfinite(ellsize)))
    warning('Non-finite projection size !!');
    elliptic_img = [];

    return;
  end

  done = false;
  while (~done)
    try
      elliptic_img = zeros(ellsize);
      done = true;
    catch ME
      warning(['The resolution of the elliptical projection is reduced as it contains more elements than the maximum available memory (' num2str(prod(ellsize)) ')']);
      ellsize = ceil(ellsize / 10);
      ellsize(ellsize < 1) = 1;
    end
  end

  if (isempty(axes_length) | any(axes_length == 0))
    return;
  end

  row = [0:ellsize(2)-1].';
  col = ones(size(row));

  for i=1:ellsize(1)
    
    [x, y] = pixels2elliptic(col * i, row, ellsize, axes_length, safety, type);
    [x, y] = elliptic2carth(x, y, center, axes_length, orient, type);

    elliptic_img(i,:) = bilinear(img,x,y);
  end

  return;
end

function [center, axes_length, orient, safety, ellsize, circular, type] = parse_inputs(varargin)

  center = [];
  axes_length = [];
  orient = [];
  safety = 1.2;
  ellsize = [];
  circular = false;
  type = 'hyperbolic';

  for i=1:length(varargin)
    var_type = get_type(varargin{i});
    switch var_type
      case 'num'
        if (isempty(center))
          center = varargin{i};
        elseif (isempty(axes_length))
          axes_length = varargin{i};
        elseif (isempty(orient))
          orient = varargin{i};
        elseif (numel(varargin{i}) == 1)
          safety = varargin{i};
        else
          ellsize = varargin{i};
        end
      case 'char'
        type = varargin{i};
      case 'bool'
        circular = varargin{i};
    end
  end

  if (isempty(axes_length))
    return;
  end

  if (length(circular) == 1)
    circular = [circular circular];
  end

  if(isempty(ellsize))
    perif = ceil(ellipse_circum(axes_length * safety)-1);
    width = ceil(axes_length(1) * safety) + 1;
    ellsize = [perif, width];

    nums = prod(ellsize);
    if (nums > 1e7)
      ellsize = floor(ellsize / (nums / 1e7));
    end
  end

  ellsize(ellsize < 1) = 1;

  return;
end
