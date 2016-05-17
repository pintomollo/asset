function [x,y] = draw_ellipse(center,axes_length,orientation,npts)

  if(nargin<4)
    npts=128;
    color = 'b';
    if (nargin == 1)
      orientation = center(5);
      axes_length = center(3:4);
      center = center(1:2);
    end
  elseif(ischar(npts))
    color = npts;
    npts = 128;
  end

  t = [0:2*pi/npts:2*pi];
  r = ones(1,npts+1);

  [ex,ey] = elliptic2carth([t;r], center, axes_length, orientation, 'radial');

  if(nargout==0)
    plot(ex,ey,color);
  elseif (nargout==1)
    x = [ex ey];
  else
    x = ex;
    y = ey;
  end
end
