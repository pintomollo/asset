function [carth_img] = carthesian_coordinate(ell_img, varargin)
  [center, axes_length, orient, safety, imgsize, circular, type] = parse_inputs(varargin{:});

  ell_size = size(ell_img);
  carth_img = zeros(imgsize);

  col = [1:imgsize(2)].';
  row = ones(size(col));

  for i=1:imgsize(1)

    [theta,rads] = carth2elliptic(col, row * i, center, axes_length, orient, type);
    [theta,rads] = elliptic2pixels(theta, rads, ell_size, axes_length, safety, type);
    
    carth_img(i,:) = bilinear(ell_img, rads, theta, double(circular));
  end

  return;
end

function [center, axes_length, orient, safety, imgsize, circular, type] = parse_inputs(varargin)

  center = [];
  axes_length = [];
  orient = [];
  safety = 1.2;
  imgsize = [];
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
          imgsize = varargin{i};
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

  return;
end
