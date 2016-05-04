function [img, minval, maxval] = imnorm(img, minval, maxval, mode, minrange, maxrange)
% IMNORM normalizes the pixel value of an image.
%
%   [IMG] = IMNORM(ORIG) normalizes the pixel value of ORIG such that it ranges
%   between 0 and 1 in IMG. A linear scaling between the min and max values of ORIG
%   is applied. Non-finite values are replaced by NaN.
%
%   [...] = IMNORM(STACK) normalizes the entire stack as a whole.
%
%   [IMG, MIN, MAX] = IMNORM(...) returns in addition the MIN and MAX values found.
%
%   [...] = IMNORM(ORIG, MIN, MAX) normalizes the values between MIN and MAX in ORIG
%   such that it ranges between 0 and 1 in IMG. Values outside [MIN, MAX] are clipped.
%
%   [...] = IMNORM(ORIG, MIN, MAX, MODE) performs the normalization according to the
%   defined MODE: either 'column', 'row' or 'slice' -wise. Provide empty values for
%   MIN and MAX to utilize the extrema values in ORIG.
%
%   [...] = IMNORM(ORIG, MIN, MAX, MODE, LBOUND, UBOUND) performs the normalization
%   such that the pixel value in IMG ranges between LBOUND and UBOUND. Provide an
%   empty string for MODE to ignore this parameter.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 15.05.2014

  % Input checking and default values
  if (nargin == 1)
    minval = [];
    maxval = [];
    mode = '';
  elseif (nargin < 4)
    mode = '';
  end

  % In case we have no data
  if (nargin == 0 || isempty(img))
    return;
  end

  % The size of the image
  [h, w, p] = size(img);

  % Get the type of img as we need to convert it to double for the division
  class_type = class(img);
  is_double = (class_type(1) == 'd');

  % In that case, convert it
  if (~is_double)
    img = double(img);
  end

  % More input checking
  if (nargin < 6)
    if (is_double)
      minrange = 0;
      maxrange = 1;
    else
      minrange = intmin(class_type);
      maxrange = intmax(class_type);
    end
  end

  % And where there are non-finite elements, we replace them with NaN which do not
  % conflict with min() and max()
  img(~isfinite(img)) = NaN;

  % Convert to lower case, just in case !
  mode = lower(mode);

  % In case we do not have a minimal value, we need to compute it
  if (isempty(minval))

    % The simplest case, we get the overall minimum
    if (isempty(mode))
      minval = min(img(:));
    else

      % Otherwise, we get it along the corresponding dimension. In addition, we
      % need to reconstruct a matrix of the proper size for the normalization.
      switch mode(1)
        case 'r'
          minval = repmat(min(img,[],2), [1, w, 1]);
        case 'c'
          minval = repmat(min(img,[],1), [h, 1, 1]);
        case 's'
          minval = repmat(min(img,[],3), [1, 1, p]);
        case 'a'
          minval = min(img(:));

        % Here we have a problem, so notify it and ignore the provided mode
        otherwise
          warning('CAST:imnorm', ['Normalization mode ' mode ' is unknown. Ignoring it.']);
          minval = min(img(:));

          % And remove the weird mode
          mode = 'a';
      end
    end
  end

  % Here we have exactly the same, except for the max value
  if (isempty(maxval))
    if (isempty(mode))
      maxval = max(img(:));
    else
      switch mode(1)
        case 'r'
          maxval = repmat(max(img,[],2), [1, w, 1]);
        case 'c'
          maxval = repmat(max(img,[],1), [h, 1, 1]);
        case 's'
          maxval = repmat(max(img,[],3), [1, 1, p]);
        case 'a'
          maxval = max(img(:));

        % Here we have a problem, so notify it and ignor the provided mode
        otherwise
          warning('CAST:imnorm', ['Normalization mode ' mode ' is unknown. Ignoring it.']);
          maxval = max(img(:));
      end
    end
  end

  % Cast all variables to double just in case
  minval = double(minval);
  maxval = double(maxval);
  minrange = double(minrange);
  maxrange = double(maxrange);

  % Now for the normalization per se. Because we kept matrices of the proper size,
  % we can compute it element-wise.
  img = (img - (minval - minrange)) .* ((maxrange - minrange) ./ (maxval - minval));

  % Clip the values outside of the defined range
  img(img < minrange) = minrange;
  img(img > maxrange) = maxrange;

  % And convert back to the original type
  if (~is_double)
    img = cast(img, class_type);
  end

  return;
end
