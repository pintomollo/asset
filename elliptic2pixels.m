function [ptsi, ptsj] = elliptic2pixels(theta, rads, imgsize, safety)

  inv = false;
  if (nargin < 4)
    [imgsize, safety] = deal(rads, imgsize);

    if(any(size(theta)==1))
      rads = theta(:);
      theta = 2 * pi * ([1:length(theta)].' - 1) / length(theta);
    else
      if (size(theta, 1) <= 2)
        theta = theta.';
        inv = true;
      end
      rads = theta(:,2);
      theta = theta(:,1);
    end
  end
  
  ptsi = (theta * imgsize(1) / (2 * pi)) + 1;
  ptsj = (rads * (imgsize(2) - 1) / safety) + 1;

  if (nargout == 1)
    ptsi = [ptsi, ptsj];
    ptsj = [];

    if (inv)
      ptsi = ptsi.';
    end
  elseif (inv)
      ptsi = ptsi.';
      ptsj = ptsj.';
  end

  return;
end
