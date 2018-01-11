function [polar_img] = polar_coordinate(img, center, orient)

  [h,w] = size(img);
  if(nargin<2)
    orient = 0;
    center = round([h,w] / 2);
  elseif(nargin<3)
    orient  = 0;
  else
    orient = (orient*180)/pi;
  end

  polar_img = zeros(2*h,w);
  diags = (sqrt(w^2+h^2)/2)/w;

  Dtheta = (1*pi)/h;
  for i=1:2*h
    theta = Dtheta * (i-1) + orient;
    for j=1:w
      radius = j * diags; 

      x = radius * cos(theta) + center(1);
      y = radius * sin(theta) + center(2);

      xf = floor(x);
      yf = floor(y);
      xc = ceil(x);
      yc = ceil(y);

      if(xf>w || yf>h || xc<1 || yc<1 || xc>w || xf<1 || yc>h || yf<1)

        polar_img(i,j) = 0;

      else
        if(xc==xf && yc==yf)
          polar_img(i,j) = img(y,x);
        elseif(xc==xf)
          polar_img(i,j) = img(yf,x)*(yc-y) + img(yc,x)*(y-yf);
        elseif(yc==yf)
          polar_img(i,j) = img(y,xf)*(xc-x) + img(y,xc)*(x-xf);
        else
          polar_img(i,j) = img(yf,xf)*(yc-y)*(xc-x) + img(yf,xc)*(yc-y)*(x-xf) + img(yc,xf)*(y-yf)*(xc-x) + img(yc,xc)*(y-yf)*(x-xf);
        end
      end
    end
  end

  return;
end
