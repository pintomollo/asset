function [img, params] = scaled_cast(img, params, type)
% ALL2UINT16 converts any type of array to uint16, rescaling it to fit the new range
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

  if (nargin == 1)
    params = [];
    type = 'uint16';
  elseif (nargin == 2)
    if (ischar(params))
      type = params;
      params = [];
    else
      type = 'uint16';
    end
  elseif (ischar(params))
    tmp = params;
    params = type;
    type = tmp;
  end

  % Check if we need to measure the parameters if IMG
  if (isempty(params))
    params = get_infos(img, type);
  end

  % If the scaling is 1, then we do not need to do anything
  if (isempty(params.target_type))
    return;
  end

  % Convert the image to double to work correctly
  img = double(img);

  % Rescale it
  img = ((img - params.offset(1)) * params.scaling) + params.offset(2);

  % Convert it to UINT16
  img = cast(img, params.target_type);

  return;
end

% For the rescaling to work, we need ot know the current range of values and
% the available one in UINT16. In addition, we store whether the type is signed or not.
function infos = get_infos(img, target)

  % Get the structure
  infos = get_struct('image_infos');

  % Get the type of numerical variable we have
  type = class(img);

  % Special case if both types are identical
  if (strncmp(target, type, length(target)))
    infos.target_type = '';
    return;
  end

  % Switch over it to get the min, max and sign information
  switch type
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
        warning('Matlab:all2uint16','Real value arrays will be rescaled');
      case {'uint8', 'uint16', 'uint32', 'uint64'}
        minval = intmin(type);
        maxval = intmax(type);
        is_signed = false;
      otherwise
        error('Matlab:all2uint16',['Unkown image data type:' class(img)]);
        return;
  end

  % Convert the parameters as Matlab do not combine different types
  minval = double(minval);
  maxval = double(maxval);

  % Switch over it to get the min, max and sign information
  switch target
      case {'int8', 'int16', 'int32', 'int64'}
        mintrg = intmin(target);
        maxtrg = intmax(target);
      case {'single', 'double'}
        mintrg = minval;
        maxtrg = maxval;
      case {'uint8', 'uint16', 'uint32', 'uint64'}
        mintrg = intmin(target);
        maxtrg = intmax(target);
      otherwise
        error('Matlab:all2uint16',['Unkown target data type:' target]);
        return;
  end

  % Convert the parameters as Matlab do not combine different types
  mintrg = double(mintrg);
  maxtrg = double(maxtrg);

  % Compute the two values
  offset = [minval mintrg];
  scaling = (maxtrg - mintrg) / (maxval - minval);

  % Store everything
  infos.is_signed = is_signed;
  infos.offset = offset;
  infos.scaling = scaling;
  infos.target_type = target;

  return;
end
