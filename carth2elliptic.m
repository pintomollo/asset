function [ellpts, radius] = carth2elliptic(varargin)

  [ptsx, ptsy, center, axes_length, orient, align, type] = parse_inputs(varargin{:});

  if (isempty(ptsx))
    ellpts = [];
    if (nargout == 2)
      radius = [];
    end

    return;
  end

  if (size(ptsy, 2) > 1)
    results = NaN(size(ptsy));
    [pts_o, results(:,1)] = carth2elliptic(ptsx, ptsy(:,1), center, axes_length, orient, align, type);

    for i=2:2:(size(ptsy, 2) - 1)
      [results(:, i:i+1)] = carth2elliptic(ptsy(:,i:i+1), [0; 0], axes_length, orient, align, type);
    end

    if (nargout == 1)
      ellpts = [pts_o results];
    else
      ellpts = [pts_o results(:,[2:2:end])];
      radius = results(:,[1:2:end]);
    end

    return;
  end

  ellpts = zeros(length(ptsx),2);

  x = ptsx - center(1);
  y = -(ptsy - center(2));

  corient = cos(orient); 
  sorient = sin(orient); 

  tmp = (x*corient + y*sorient);
  y = (y*corient - x*sorient);
  x = tmp;

  switch type
    case 'hyperbolic'
      focus = sqrt(-diff(axes_length.^2));

      ga = sqrt((x + focus).^2 + y.^2);
      gb = sqrt((x - focus).^2 + y.^2);

      signs = sign(y);
      signs(signs==0) = 1;

      ellpts(:,1) = signs .* acos((ga - gb) / (2*focus));
      ellpts(:,2) = acosh((ga + gb) / (2*focus));

    case 'radial'
      ratio = axes_length(1) / axes_length(2);

      long_axe = x/axes_length(1);
      short_axe = y/axes_length(2);
      
      ellpts(:,1) = atan2(short_axe, long_axe);
      ellpts(:,2) = hypot(long_axe, short_axe);
    otherwise
      warning(['Unkown projection type: "' type '", using "hyperbolic" instead.']);
      [ellpts] = carth2elliptic(ptsx, ptsy, center, axes_length, orient, align, 'hyperbolic');
  end

  ellpts = real(ellpts);
  negs = ellpts(:,1)<0;
  ellpts(negs,1) = ellpts(negs,1) + 2*pi;

  % Correct a rounding error
  if (abs(ellpts(1,1) - 2*pi) < 1e-6)
    ellpts(1,1) = 0;
  end

  if (align)
    indx = find(ellpts(:,1)==min(ellpts(:,1)));

    if (all(ellpts(1,:)) == ellpts(end,:))
      ellpts = [ellpts(indx:end-1,:); ellpts(1:indx-1,:); ellpts(indx,:)];
    else
      ellpts = [ellpts(indx:end,:); ellpts(1:indx-1,:)];
    end
  end

  if(nargout==2)
    radius = ellpts(:,2);
    ellpts(:,2) = [];
  end

  return;
end

function [ptsx, ptsy, center, axes_length, orient, align, type] = parse_inputs(varargin)

  ptsx = [];
  ptsy = [];
  center = [];
  axes_length = [];
  orient = [];
  align = false;
  type = 'hyperbolic';

  for i=1:length(varargin)
    var_type = get_type(varargin{i});
    switch var_type
      case 'num'
        if (isempty(ptsx))
          ptsx = varargin{i};
        elseif (all(size(ptsx) == size(varargin{i})))
          ptsy = varargin{i};
        elseif (isempty(center) & numel(varargin{i})==2)
          center = varargin{i};
        elseif (isempty(axes_length) & numel(varargin{i})==2)
          axes_length = varargin{i};
        elseif (isempty(orient) & numel(varargin{i})==1)
          orient = varargin{i};
        end
      case 'char'
        type = varargin{i};
      case 'bool'
        align = varargin{i};
      case 'none'
        if (isempty(ptsx))
          ptsx = NaN(1,2);
        end
    end
  end

  if (isempty(ptsx) | all(isnan(ptsx)))
    ptsx = [];
    return;
  elseif (size(ptsx,2) > 4)
    ptsx = ptsx.';
  end

  if (isempty(ptsy))
    ptsy = ptsx(:, 2:end);
    ptsx = ptsx(:,1);
  end

  if (isempty(center))
    [center, axes_length, orient] = fit_ellipse(ptsx, ptsy(:, 1));
  end

  return;
end
