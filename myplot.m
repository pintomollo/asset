function handles = myplot(pts, varargin)

  [h, polar, args] = parse_inputs(varargin{:});

  if (iscell(pts))
    pts = pts{1};
  end

  [m, n] = size(pts);
  if (m <= 4 & n > 4)
    pts = pts.';

    [m, n] = size(pts);
  end

  if (n == 1)
    if (length(args) > 0 & isnumeric(args{1}))
      tmp = args{1};
      args(1) = [];

      tmp = reshape(tmp, size(pts));
      pts = [pts tmp];

      [m, n] = size(pts);
    else
      pts = [pts [1:length(pts)].'];

      [m, n] = size(pts);
    end
  end

  for i=1:length(args)
    if (ischar(args{i}) & length(args{i}) == 1)
      switch args{i}
        case 'r'
          args = [args(1:i-1) {'Color'} {[1 0 0]} args(i+1:end)];
        case 'g'
          args = [args(1:i-1) {'Color'} {[0 1 0]} args(i+1:end)];
        case 'b'
          args = [args(1:i-1) {'Color'} {[0 0 1]} args(i+1:end)];
        case 'y'
          args = [args(1:i-1) {'Color'} {[1 1 0]} args(i+1:end)];
        case 'm'
          args = [args(1:i-1) {'Color'} {[1 0 1]} args(i+1:end)];
        case 'c'
          args = [args(1:i-1) {'Color'} {[0 1 1]} args(i+1:end)];
        case 'w'
          args = [args(1:i-1) {'Color'} {[1 1 1]} args(i+1:end)];
        case 'k'
          args = [args(1:i-1) {'Color'} {[0 0 0]} args(i+1:end)];
      end
    end
  end

  if (isempty(h))
    h = line(pts(:,1), pts(:,2), args{:});
  else
    set(h(1), 'XData', pts(:,1), 'YData', pts(:,2), args{:});
  end

  if (n > 2)

    if (n == 3)
      pts = pts(:,[1:end end]);
    end

    indx = 0;
    for i=1:2:length(args)
      if (strcmp(args{i}, 'Color'))
        indx = i;
      end
    end

    if (indx > 0)
      args(indx:indx+1) = [];
    end 

    color = get(h(1), 'Color');
    color = [ones(size(color)); color];
    color = mean(color, 1);

    if (numel(h) == 1)
      h(end+1) = surface('XData', [pts(:,1) - pts(:,3) pts(:,1) + pts(:,3)], ...
                         'YData', [pts(:,2) - pts(:,4) pts(:,2) + pts(:,4)], ...
                         'ZData', zeros(m, 2), ...
                         args{:}, ...
                         'FaceColor', color, ...
                         'LineStyle', 'none', ...
                         'Marker', 'none');
    else
      set(h(2),'XData', [pts(:,1) - pts(:,3) pts(:,1) + pts(:,3)], ...
               'YData', [pts(:,2) - pts(:,4) pts(:,2) + pts(:,4)], ...
               'ZData', zeros(m, 2), ...
               args{:}, ...
               'FaceColor', color, ...
               'LineStyle', 'none', ...
               'Marker', 'none');
    end
    
    uistack(h(end), 'down');
  else
    h = h(1);
  end

  if (nargout > 0)
    handles = h;
  end

  return;
end

function [handles, polar, args] = parse_inputs(varargin)
  
  handles = [];
  polar = false;
  args = {};

  if (nargin > 0)
    for i=1:length(varargin)
      type = get_type(varargin{i});
      switch(type)
        case 'bool'
          polar = varargin{i};
        case 'num'
          if (all(ishandle(varargin{i})))
            handles = varargin{i};
          end
        case 'char'
          args = [args varargin(i:end)];
          return;

        otherwise
          args = [args varargin(i)];
      end
    end
  end

  return;
end
