function [theta,rads] = pixels2elliptic(ptsi, ptsj, imgsize, safety)

  if (nargin < 4)
    [imgsize, safety] = deal(ptsj, imgsize);

    if(size(ptsi,2)==1)
      ptsj = ptsi;
      ptsi = [1:length(ptsi)].';
    else
      ptsj = ptsi(:,2);
      ptsi = ptsi(:,1);
    end
  end

  theta = 2 * pi * (ptsi - 1) / imgsize(1);
  rads = safety * (ptsj - 1) / (imgsize(2) - 1);

  if (nargout == 1)
    theta = [theta, rads];
    rads = [];
  end

  return;
end
