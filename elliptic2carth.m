function [carth_pts, carth_y] = elliptic2carth(varargin)

  [ptso, ptsr, center, axes_length, orient, type] = parse_inputs(varargin{:});
  
  if (numel(ptso)==0)
    carth_pts = ptso;

    if(nargout==2)
      carth_y = ptso;
    end

    return;
  end

  if (size(ptsr, 2) > 1)
    results = NaN(size(ptsr));
    [pts_x, results(:,1)] = elliptic2carth(ptso, ptsr(:,1), center, axes_length, orient, type);

    for i=2:2:(size(ptsr, 2) - 1)
      [results(:, i:i+1)] = elliptic2carth(ptsr(:,i:i+1), [0; 0], axes_length, orient, type);
    end

    if (nargout == 1)
      carth_pts = [pts_x results];
    else
      carth_pts = [pts_x results(:,[2:2:end])];
      carth_y = results(:,[1:2:end]);
    end

    return;
  end

  orient = orient; 
  corient = cos(orient);
  sorient = sin(orient);

  O = ptso;
  r = ptsr;

  carth_pts = zeros(length(ptso),2);

  switch type 
    case 'hyperbolic'

      focus = sqrt(-diff(axes_length.^2));
      carth_pts = zeros(length(ptso),2);
      carth_pts(:,1) = focus*cosh(r).*cos(O)*corient - ...
                       focus*sinh(r).*sin(O)*sorient + center(1);
      carth_pts(:,2) = -(focus*cosh(r).*cos(O)*sorient + ...
                         focus*sinh(r).*sin(O)*corient) + center(2);
    case 'radial'

      carth_pts(:,1) = axes_length(1)*r.*cos(O)*corient - ...
                       axes_length(2)*r.*sin(O)*sorient + center(1);
      carth_pts(:,2) = -(axes_length(1)*r.*cos(O)*sorient + ...
                         axes_length(2)*r.*sin(O)*corient) + center(2);

    otherwise
      warning(['Unkown projection type: "' type '", using "hyperbolic" instead.']);
      [carth_pts] = elliptic2carth(O, r, center, axes_length, orient, 'hyperbolic');
  end

  if(nargout==2)
    carth_y = carth_pts(:,2);
    carth_pts(:,2) = [];
  end

  return;
end

function [ptso, ptsr, center, axes_length, orient, type] = parse_inputs(varargin)

  ptso = [];
  ptsr = [];
  center = [];
  axes_length = [];
  orient = [];
  type = 'hyperbolic';

  for i=1:length(varargin)
    var_type = get_type(varargin{i});
    switch var_type
      case 'num'
        if (isempty(ptso))
          ptso = varargin{i};
        elseif (isempty(ptsr) & ndims(ptso) == ndims(varargin{i}) & all(size(ptso) == size(varargin{i})))
          ptsr = varargin{i};
        elseif (isempty(center) & numel(varargin{i})==2)
          center = varargin{i};
        elseif (isempty(axes_length) & numel(varargin{i})==2)
          axes_length = varargin{i};
        elseif (isempty(orient) & numel(varargin{i})==1)
          orient = varargin{i};
        end
      case 'char'
        type = varargin{i};
      case 'none'
        if (isempty(ptso))
          ptso = NaN(1,2);
        end
    end
  end

  if (isempty(ptso) | all(isnan(ptso)))
    ptso = [];
    return;
  elseif (size(ptso,2) > 4)
    ptso = ptso.';
  end

  if (isempty(ptsr))
    ptsr = ptso(:, 2:end);
    ptso = ptso(:,1);
  end

  if (isempty(center))
    ptsr = [];
    ptso = [];
  end

  return;
end
