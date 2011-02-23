function [carth_img] = carthesian_coordinate(ell_img, center, axes_length, orient, safety, imgsize, circular)

  if (nargin < 7)
    circular = false;
  end

  [perif, width] = size(ell_img);
  carth_img = zeros(imgsize);

  col = [1:imgsize(2)].';
  row = ones(size(col));

  for i=1:imgsize(1)

    [theta,rads] = carth2elliptic(col, row * i, center, axes_length, orient);
    [theta,rads] = elliptic2pixels(theta, rads, [perif, width], safety);
    
    carth_img(i,:) = bilinear_mex(ell_img, rads, theta, circular);
    %carth_img(i,:) = bilinear(ell_img,rads,theta, circular);
  end

  return;
end
