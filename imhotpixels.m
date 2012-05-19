function img = imhotpixels(orig_img, varargin)
% IMHOTPIXELS removes the so-called "hot pixels" in an image. These pixels are defined
% as having a value larger/smaller than MEAN(pixels) +/- THRESH*STD(pixels). The 
% pixels detected as such are then replaced by another value which is computed using 
% the neighboring values.
%
%   [IMG] = IMHOTPIXELS(IMG) filters IMG using the default parameters i.e. THRESH = 15
%   and the replacement method is the MEDIAN function.
%
%   [IMG] =  IMHOTPIXELS(IMG, THRESH) specifies the value of THRESH.
%
%   [IMG] =  IMHOTPIXELS(IMG, METHOD, PARAMS) specifies the type of method used to 
%   compute the new value of the hot pixels. METHOD should be a string which value
%   correspond to one of the fspecial filter type (help fspecial). PARAMS contains 
%   the parameters of the fspecial filter as a cell array if more than one parameter
%   is required. If none is required, PARAMS should be an empty cell.
%
%   [IMG] = IMHOTPIXELS(IMG, 'custom', FILTER) uses a custom 2D-filter to compute the
%   new value of the hot pixels, which is defined by FILTER.
%
%   [IMG] = IMOTPIXELS(IMG, THRESH, METHOD, PARAMS) specifies both.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 30.06.2011

  % Sort out the different types of input and assign default values
  [thresh, method, params, filter] = parse_input(varargin{:});

  % Get the pixels as a vector
  pixels = double(orig_img(:));
  % Compute their mean and standard deviation
  [mean_value, stddev] = mymean(pixels);

  % Get the threshold value
  thresh = thresh * stddev;

  % Identify the hot pixels
  bad_pixels = ((orig_img > (mean_value + thresh)) | (orig_img < mean_value - thresh));

  % Create the output image
  img = orig_img;

  % If anything has to be done
  if (any(bad_pixels(:)))
  
    % Apply the required strategy to compute the new value of the hot pixels
    switch (method)

      % Plain median filtering
      case 'median'
        filt_img = median_mex(double(orig_img));

      % Using the provided custom filter
      case 'custom'
        filt_img = imfilter(double(orig_img), filter, 'symmetric');

      % Using the specified fspecial filter
      otherwise
        if (isempty(params))
          filt_img = imfilter(double(orig_img), fspecial(method), 'symmetric');
        else
          filt_img = imfilter(double(orig_img), fspecial(method, params{:}), 'symmetric');
        end
    end

    % Replace the hot pixels
    img(bad_pixels) = filt_img(bad_pixels);
  end

  return;
end

% Sorts out the different inputs
function [thresh, method, params, filter] = parse_input(varargin)

  % The default values
  thresh = 15;
  method = 'median';
  params = [];
  filter = [];

  % Loop over all the inputs
  for i = 1:length(varargin)

    % Assign them according to their type
    switch (get_type(varargin{i}))

      % Must be the method type
      case 'char'
        method = varargin{i};

      % Parameters for the fspecial filter
      case 'cell'
        params = varargin{i};

      % Several possibilities
      case 'num'

        % That many values is the cutom filter
        if (numel(varargin{i}) > 2)
          filter = varargin{i};
          method = 'custom';

        % Single value with the median method which has no parameter is the threshold
        elseif (numel(varargin{i}) == 1 & strncmp(method, 'median', 6))
          thresh = varargin{i};

        % Otherwise it's the parameters
        else
          params = varargin{i};
        end
    end
  end

  % We use only cells for the parameters
  if (~iscell(params))
    params = {params};
  end
  
  return;
end
