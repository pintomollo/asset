function [img, params] = all2uint16(img, params)
% ALL2UINT16 converts any type of array to uint16, rescaling it ot fit the new range
% of values.
%
%   M = ALL2UINT16(M) converts M to UINT16.
%
%   [M, PARAMS] = ALL2UINT16(M) returns some parameters related to the original M.
%
%   [...] = ALL2UINT16(M, PARAMS) provides the parameters related to M, this is
%   useful to speed up batch processing (e.g. for stack of images). If PARAMS is
%   empty, it will be measured.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 20.06.2011

  % Check if we need to measure the parameters if IMG
  if (nargin == 1 | isempty(params))
    params = get_infos(img);
  end

  % If the scaling is 1, then we do not need to do anything
  if (params.scaling == 1)
    return;
  end

  % Convert the image to double to work correctly
  img = double(img);

  % Rescale it
  img = (img + params.offset) * params.scaling;

  % Convert it to UINT16
  img = uint16(img);

  return;
end

% For the rescaling to work, we need ot know the current range of values and
% the available one in UINT16. In addition, we store whether the type is signed or not.
function infos = get_infos(img)

  % Get the type of numerical variable we have
  type = class(img);

  % Switch over it to get the min, max and sign information
  switch class(img)
      case {'int8', 'int16', 'int32', 'int64'}
        minval = intmin(type);
        maxval = intmax(type);
        is_signed = true;
      case {'single', 'double'}
        minval = min(img(:));
        maxval = max(img(:));
        is_signed = true;

        % Warn the user as there is no "range" for double so we will use the
        % range of values
        warning('Real value arrays will be rescaled');
      case {'uint8', 'uint16', 'uint32', 'uint64'}
        minval = intmin(type);
        maxval = intmax(type);
        is_signed = false;
  end

  % Get the futur range of values
  max16 = double(intmax('uint16'));

  % Convert the parameters as Matlab do not combine different types
  minval = double(minval);
  maxval = double(maxval);
  
  % Compute the two values
  offset = -minval;
  scaling = max16 / (maxval - minval);

  % Create the structure
  infos = struct('is_signed', is_signed, ...
                 'offset', offset, ...
                 'scaling', scaling);

  return;
end
