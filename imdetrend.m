function [img, p] = imdetrend(orig_img, npts)

  if (nargin == 1)
    npts = 16;
  end

  [h,w] = size(orig_img);
  x = floor(npts/2):floor(w/npts):w;
  y = floor(npts/2):floor(h/npts):h;

  img = imopen(orig_img, strel('disk',10));
  img = gaussian_mex(img, 10);
  %img = imfilter(img, fspecial('gaussian', 20, 10), 'symmetric');

  trend = img(y,x);

  [X, Y] = meshgrid(x,y);

  indx = (X < w/4 | X > 3*w/4 | Y < h/4 | Y > 3*h/4);

  p = [X(indx).^2 X(indx).*Y(indx) Y(indx).^2 X(indx) Y(indx) ones(numel(trend(indx)),1)] \ trend(indx);

  [X, Y] = meshgrid(1:w,1:h);
  bkg = p(1) * X.^2 + p(2) * Y.*X + Y.^2 * p(3) + X * p(4) + Y * p(5) + p(6);

  img = imsubtract(orig_img, bkg);
  %img = imnorm(img);

  return;
end
