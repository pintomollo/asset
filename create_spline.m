function [pp] = create_spline(pts, varargin)
%CSAPE Cubic spline interpolation with various end conditions.

  try
  pp = get_struct('spline',1);
  [corr_orient, pos] = parse_input(varargin{:});

  if (isempty(pts))
    return;
  end

  if (size(pts,2) <= 4)
    pts = pts.';
  end
  pos = pos(:).';

  if (~isempty(pos) && size(pos,2) ~= size(pts,2))
    return;  
  end

  indxs_pts = ~any((isnan(pts) | isinf(pts)),1);
  if (~isempty(pos))
    indxs_pts = indxs_pts & ~any((isnan(pos) | isinf(pos)),1);
    pos = pos(indxs_pts);
  end
  pts = pts(:,indxs_pts);

  if (size(pts,2) < 2)
    % Only one point, so create a "circular" spline around it
    if (numel(pts) == 0)
      return;
    end

    pts = [pts pts pts];
    pos = [0 pi 2*pi];
  end

  if (isempty(pos))
    if (any(pts(:,1) ~= pts(:,end)))
      pts = pts(:,[1:end 1]);
    end

    if (corr_orient & ~ispolycw(pts(1,:),pts(2,:)))
      pts = pts(:,[end:-1:1]);
    end

    dt = sum((diff(pts.').^2).');

    indxs = logical([dt~=0 1]); %& ~(any((isnan(pts) | isinf(pts)),1));
    pts = pts(:,indxs);
    dt = dt(:,indxs(1:end-1));
    t = cumsum([0,dt.^(1/4)]);
  else

    if (corr_orient)
      if (size(pts, 1) < 2)
        [pts, pos] = poly2cw(pts, pos);
      elseif (~ispolycw(pts(1,:),pts(2,:)))
        pts = pts(:,[end:-1:1]);
        pos = pos(:,[end:-1:1]);
      end
    end

    dt = diff(pos);

    indxs = logical([dt~=0 1]); %& ~(any((isnan(pts) | isinf(pts)),1));

    pts = pts(:,indxs);
    t = pos(:,indxs);
  end

  if (isempty(pts))
    return;
  end

  if (length(t) < 2)
    % Only one point, so create a "circular" spline around it
    if (numel(pts) == 0)
      return;
    end

    pts = [pts pts pts];
    t = [0 pi 2*pi];
  end

    pp = create_really(pts, t);
  catch ME
    pts
    keyboard
  end

  return;
end

function [pp] = create_really(pts, t)
 
  %x = t;
  %y = pts;

  [xi,yi,sizeval,endvals] = chckxywp(t,pts,0);
  [n,yd] = size(yi);
  dd = ones(1,yd);
  dx = diff(xi);
  divdif = diff(yi)./dx(:,dd);

% set up the linear system for solving for the slopes at XI.
%dx = diff(xi); divdif = diff(yi)./dx(:,dd);

  c = spdiags([ [dx(2:n-1,1);0;0] ...
              2*[0;dx(2:n-1,1)+dx(1:n-2,1);0] ...
                [0;0;dx(1:n-2,1)] ], [-1 0 1], n, n);
  b = zeros(n,yd);
  b(2:n-1,:)=3*(dx(2:n-1,dd).*divdif(1:n-2,:)+dx(1:n-2,dd).*divdif(2:n-1,:));
  c(1,1)=1; c(1,n)=-1;
  c(n,1:2)=dx(n-1)*[2 1];
  c(n,n-1:n)= c(n,n-1:n)+dx(1)*[1 2];
  b(n,:) = 3*(dx(n-1)*divdif(1,:) + dx(1)*divdif(n-1,:));

  % solve for the slopes ..  (protect current spparms setting)
  %mmdflag = spparms('autommd');
  %spparms('autommd',0); % suppress pivoting
  mmdflag = sparsfun('slashset');
  sparsfun('slashset', [0 mmdflag(2:end)]);
  s=c\b;
  sparsfun('slashset', mmdflag);
  %spparms('autommd',mmdflag);

  %                          .. and convert to ppform
  c4 = (s(1:n-1,:)+s(2:n,:)-2*divdif(1:n-1,:))./dx(:,dd);
  c3 = (divdif(1:n-1,:)-s(1:n-1,:))./dx(:,dd) - c4;
  pp = ppmak(xi.', ...
   reshape([(c4./dx(:,dd)).' c3.' s(1:n-1,:).' yi(1:n-1,:).'],(n-1)*yd,4),yd);

  return;
end

function [corr_orient, pos] = parse_input(varargin)

  corr_orient = true;
  pos = [];

  if (nargin > 0)
    for i=1:length(varargin)
      type = get_type(varargin{i});
      switch type
        case 'bool'
          corr_orient = varargin{i};
        case 'num'
          pos = varargin{i};
      end
    end
  end

  return;
end
