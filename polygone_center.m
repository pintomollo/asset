function [center, area] = polygone_center(pts)
  
  if (size(pts,1)==2)
    pts = pts.';
  end

  if (~all(pts(1,:)==pts(end,:)))
    pts = pts([1:end 1],:);
  end

  npts = size(pts,1);
  indx = 1:npts-1;
  next = 2:npts;

  tmp = pts(indx,1).*pts(next,2) - pts(next,1).*pts(indx,2);
  area = sum(tmp) / 2;

  if (area == 0)
    center = pts(1,:);
  else
    center = (sum((pts(indx,:)+pts(next,:)) .* tmp(:,[1 1])) / (6 * area));
  end

  if (nargout<2)
    area = [];
  end

  return;
end
