function vals = bilinear(img, x, y, circular)

  if (nargin < 3)
    y = x(2,:);
    x = x(1,:);

    circular = false;
  elseif (nargin < 4)
    if (islogical(y))
      circular = y;

      y = x(2,:);
      x = x(1,:);
    else
      circular = false;
    end
  end

  if (length(circular) == 1)
    circular = [circular circular];
  end

  if (exist('bilinear_mex') == 3)
    vals = bilinear_mex(img, x, y, circular);

    return;
  end

  vals = NaN(size(x));

  [h,w] = size(img);

  xf = floor(x);
  yf = floor(y);
  xc = ceil(x);
  yc = ceil(y);

  dxf = x - xf;
  dyf = y - yf;
  dxc = xc - x; 
  dyc = yc - y;

  if (circular(1))
    xf = mod(xf,w);
    xc = mod(xc,w);

    xf(xf==0) = w;
    xc(xc==0) = w;
  end

  if (circular(2))
    yf = mod(yf,h);
    yc = mod(yc,h);

    yf(yf==0) = h;
    yc(yc==0) = h;
  end

  out = (xf>w | yf>h | xc<1 | yc<1 | xc>w | xf<1 | yc>h | yf<1);

  single = (xc==xf & yc==yf) & ~out;
  vals(single) = img(sub2ind([h,w], y(single), x(single)));
  
  vert = (xc==xf) & ~out & ~single;
  vals(vert) = img(sub2ind([h,w], yf(vert),x(vert))).*dyc(vert) + img(sub2ind([h,w], yc(vert),x(vert))).*dyf(vert);
    
  horiz = (yc==yf) & ~out & ~single;
  vals(horiz) = img(sub2ind([h,w], y(horiz),xf(horiz))).*dxc(horiz) + img(sub2ind([h,w], y(horiz),xc(horiz))).*dxf(horiz);
  
  full = ~(out | single | vert | horiz);
  vals(full) = img(sub2ind([h,w], yf(full),xf(full))).*dyc(full).*dxc(full) + ...
               img(sub2ind([h,w], yf(full),xc(full))).*dyc(full).*dxf(full) + ...
               img(sub2ind([h,w], yc(full),xf(full))).*dyf(full).*dxc(full) + ...
               img(sub2ind([h,w], yc(full),xc(full))).*dyf(full).*dxf(full);
end

