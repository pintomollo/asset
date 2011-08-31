function [theta,rads] = pixels2elliptic(varargin)

  [ptsi, ptsj, imgsize, axes_length, safety, type] = parse_inputs(varargin{:});

  theta = zeros(length(ptsi),2);

  theta(:,1) = 2 * pi * (ptsi - 1) / imgsize(1);
  theta(:,2) = safety * (ptsj - 1) / (imgsize(2) - 1);

  switch type 
    case 'radial'
    case 'hyperbolic'
      focus = sqrt(-diff(axes_length.^2));

      scaling = acosh(axes_length(1) * safety / focus);
      theta(:,2) = theta(:, 2) * (scaling / safety);

    otherwise
      warning(['Unkown projection type: "' type '", using "hyperbolic" instead.']);
      theta = pixels2elliptic(ptsi, ptsj, imgsize, axes_length, safety, 'hyperbolic');
  end

  if (nargout == 2)
    rads = theta(:,2);
    theta = theta(:,1);
  end

  return;
end

function [ptsi, ptsj, imgsize, axes_length, safety, type] = parse_inputs(varargin)

  ptsi = [];
  ptsj = [];
  imgsize = [];
  axes_length = [];
  safety = 1.2;
  type = 'hyperbolic';

  for i=1:length(varargin)
    var_type = get_type(varargin{i});
    switch var_type
      case 'num'
        if (isempty(ptsi))
          ptsi = varargin{i};
        elseif (all(size(ptsi) == size(varargin{i})))
          ptsj = varargin{i};
        elseif (isempty(imgsize) & numel(varargin{i})==2)
          imgsize = varargin{i};
        elseif (numel(varargin{i}) == 1)
          safety = varargin{i};
        else
          axes_length = varargin{i};
        end
      case 'char'
        type = varargin{i};
      case 'none'
        if (isempty(ptsi))
          ptsi = NaN(1,2);
        end
    end
  end

  if (isempty(ptsi))
    return;
  elseif (size(ptsi,2) > 4)
    ptsi = ptsi.';
  end

  if (isempty(ptsj))
    if (size(ptsi, 2) > 1)
      ptsj = ptsi(:, 2:end);
      ptsi = ptsi(:,1);
    else  
      ptsj = ptsi;
      ptsi = [1:length(ptsi)].';
    end
  end

  if (isempty(axes_length))
    type = 'radial';
  end

  return;
end
