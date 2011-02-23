function [x,y] = draw_ellipse(center,axes_length,orientation,npts)
  
  if(nargin<4)
    npts=128;
    color = 'b';
  elseif(ischar(npts))
    color = npts;
    npts = 128;
  end

  t = [0:2*pi/npts:2*pi];
  r = ones(1,npts+1);

  [x,y] = elliptic2carth([t;r], center, axes_length, orientation);
  %x=axes_length(1)*cos(t)*cos(orientation) - axes_length(2)*sin(t)*sin(orientation) + center(1);
  %y=axes_length(1)*cos(t)*sin(orientation) + axes_length(2)*sin(t)*cos(orientation) + center(2);

  if(nargout==0)
    plot(x,y,color);
    x = [];
    y = [];
  elseif (nargout==1)
    x = [x y];
    y = [];
  end
end
