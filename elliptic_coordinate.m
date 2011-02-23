function [elliptic_img] = elliptic_coordinate(img, center, axes_length, orient, safety, ellsize, circular)

  if (nargin < 5)
    safety = 1;
    ellsize = [];
    circular = false;
  elseif (nargin < 6)
    if (islogical(safety))
      circular = safety;
      safety = 1;
      ellsize = [];
    elseif (prod(size(safety)) > 1)
      ellsize = safety;
      safety = 1;
    else
      ellsize = [];
    end
  elseif (nargin < 7)
    if (~islogical(ellsize))
      circular = false;
    elseif (prod(size(safety)) > 1)
      [ellsize, circular] = deal(safety, ellsize);
      safety = 1;
    else
      circular = ellsize;
      ellsize = [];
    end
  end

  if(isempty(ellsize))
    perif = ceil(ellipse_circum(axes_length * safety)-1);
    width = ceil(axes_length(1) * safety) + 1;
  else
    perif = ellsize(1);
    width = ellsize(2);
  end

  if (perif < 1)  
    perif = 1;
  end
  if (width < 1)
    width = 1;
  end

  [h,w] = size(img);
  done = false;
  while (~done)
    try
      elliptic_img = zeros(perif,width);
      done = true;
    catch ME
      warning(['The resolution of the elliptical projection is reduced as it contains more elements than the maximum available (' num2str(perif * width) ')']);
      perif = perif / 2;
      width = width / 2;
    end
  end

  if (any(axes_length == 0))
    return;
  end

  row = [0:width-1].';
  col = ones(size(row));

  for i=1:perif
    
    [x, y] = pixels2elliptic(col * i, row, [perif, width], safety);
    [x, y] = elliptic2carth(x, y, center, axes_length, orient);

    elliptic_img(i,:) = bilinear_mex(img,x,y);
  end

  return;
end
