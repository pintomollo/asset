function [carth_img] = carthesian_coordinate(ell_img, varargin)
  [center, axes_length, orient, safety, imgsize, circular, type] = parse_inputs(varargin{:});

  ell_size = size(ell_img);
  %carth_img = zeros(imgsize);

  %col = [1:imgsize(2)].';
  %row = ones(size(col));

  col = [1:imgsize(2)];
  row = [1:imgsize(1)].';

  [I,J] = meshgrid(row, col);
  [theta,rads] = carth2elliptic(J(:), I(:), center, axes_length, orient, type);
  [theta,rads] = elliptic2pixels(theta, rads, ell_size, axes_length, safety, type);
  carth_img = bilinear(ell_img, rads, theta, double(circular));

  carth_img = reshape(carth_img, imgsize);

%  for i=1:imgsize(1)

%    [theta,rads] = carth2elliptic(col, row * i, center, axes_length, orient, type);
%    [theta,rads] = elliptic2pixels(theta, rads, ell_size, axes_length, safety, type);
    
%    carth_img(i,:) = bilinear(ell_img, rads, theta, double(circular));
%  end

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
    %var_type = get_type(varargin{i});
    if (~isempty(varargin{i}))
      var_type = class(varargin{i});
      switch var_type
        case {'double', 'single', 'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32', 'uint64'}
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
        case 'logical'
          circular = varargin{i};
      end
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
