function [ptsi, ptsj] = elliptic2pixels(varargin)

  [theta, rads, imgsize, axes_length, safety, type] = parse_inputs(varargin{:});
%function [ptsi, ptsj] = elliptic2pixels(theta, rads, imgsize, safety)

  ptsi = zeros(length(theta),2);

  ptsi(:,1) = (theta * imgsize(1) / (2 * pi)) + 1;
  switch type
    case 'radial'
      ptsi(:,2) = (rads * (imgsize(2) - 1) / safety) + 1;
    case 'hyperbolic'
      focus = sqrt(-diff(axes_length.^2));

      x = focus * cosh(rads) .* cos(theta);
      y = focus * sinh(rads) .* sin(theta);

      r = hypot(x/axes_length(1), y/axes_length(2));

      ptsi(:,2) = (r * (imgsize(2) - 1) / safety) + 1;
    otherwise
      warning(['Unkown projection type: "' type '", using "hyperbolic" instead.']);
      theta = elliptic2pixels(theta, rads, imgsize, axes_length, safety, 'hyperbolic');
  end


  if (nargout == 2)
    ptsj = ptsi(:,2);
    ptsi = ptsi(:,1);
  end

  return;
end

function [theta, rads, imgsize, axes_length, safety, type] = parse_inputs(varargin)

  theta = [];
  rads = [];
  imgsize = [];
  axes_length = [];
  safety = 1.2;
  type = 'hyperbolic';

  for i=1:length(varargin)
    var_type = get_type(varargin{i});
    switch var_type
      case 'num'
        if (isempty(theta))
          theta = varargin{i};
        elseif (numel(theta) == numel(varargin{i}))
          rads = varargin{i};
        elseif (isempty(imgsize))
          imgsize = varargin{i};
        elseif (numel(varargin{i}) == 1)
          safety = varargin{i};
        else
          axes_length = varargin{i};
        end
      case 'char'
        type = varargin{i};
    end
  end

  if (isempty(theta))
    return;
  elseif (size(theta,2) > 4)
    theta = theta.';
  end

  if (isempty(rads))
    if (size(theta, 2) > 1)
      rads = theta(:, 2:end);
      theta = theta(:,1);
    else  
      rads = theta;
      theta = 2 * pi * ([1:length(theta)].' - 1) / length(theta);
    end
  end

  if (isempty(axes_length))
    type = 'radial';
  end

  return;
end
