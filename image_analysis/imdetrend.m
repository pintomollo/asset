function [img, p] = imdetrend(orig_img, npts)
% IMDETREND removes a 2D quadratic trend in an image as produced by uneven illumination.
%
%   [IMG] = IMDETREND(IMG) removes the trend in IMG which is computed by smoothing
%   IMG and performing a 2D quadratic least-square fit onto the outter pixels of the
%   resulting image. Only the outter most fourth of the image is used for the fit as
%   the central part if often occupied by bright elements of interest.
%
%   [IMG, COEFS] = IMDETREND(...) returns in addition the 6 coefficients of the fit.
%   The coefficients are provided in the following order: [X^2 X*Y Y^2 X Y 1]
%
%   [...] = IMDETREND(IMG, NPTS) provides in addition the number of points used for
%   the fit. NPTS can be [NPTS_X NPTS_Y] to obtain an anisotropic sampling resolution.
%   The default value is 32 for both dimension.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 07.07.2011

  % Default value for the sampling
  if (nargin == 1)
    npts = 32;
  end

  % Duplicate the sampling to get one value for each dimension
  if (numel(npts) == 1)
    npts = [npts npts];
  end

  % Get the size if the image
  [h, w] = size(orig_img);
  img_size = max(h, w);

  % And deduce which pixels we'll pick to get the correct number of points
  x = floor(npts(1)/2):floor(w/npts(1)):w;
  y = floor(npts(2)/2):floor(h/npts(2)):h;

  % Construct the "grid" for the index of the pixel pairs
  [X, Y] = meshgrid(x, y);

  % Keep only the outter ones
  indx = (X < w/4 | X > 3*w/4 | Y < h/4 | Y > 3*h/4);
  X = X(indx);
  Y = Y(indx);

  % Get their index
  indx = sub2ind([h, w], Y, X);

  % Refresh the packages
  %pkg unload image
  %pkg load image

  % Filter intensively the image to keep only the trend
  img = imopen(double(orig_img), strel('disk', ceil(img_size/75), 0));
  img = gaussian_mex(img, ceil(img_size/100));

  %pkg unload image

  % Extract the pixels
  trend = img(indx);

  % Perform the least-square fit
  p = [X.^2 X.*Y Y.^2 X Y ones(size(indx))] \ trend;

  % Compute the indexes corresponding to each pixel of the image
  [X, Y] = meshgrid(1:w, 1:h);

  % And compute the corresponding trend
  bkg = p(1) * X.^2 + p(2) * Y.*X + Y.^2 * p(3) + X * p(4) + Y * p(5) + p(6);

  % Finally, substract the trend to the image
  img = double(orig_img) - bkg + mean(bkg(:));

  % Set the type of the image back to its original one
  img = cast(img, class(orig_img));

  return;
end
