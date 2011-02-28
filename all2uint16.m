function [img, params] = all2uint16(img, params)

  if (nargin == 1 | isempty(params))
    params = get_infos(img);
  end

  if (params.scaling == 1)
    return;
  end

  img = double(img);
  img = (img + params.offset) * params.scaling;
  img = uint16(img);

  return;
end

function infos = get_infos(img)

  type = class(img);

  switch class(img)
      case {'int8', 'int16', 'int32', 'int64'}
        minval = intmin(type);
        maxval = intmax(type);
        issigned = true;
      case {'single', 'double'}
        minval = min(img(:));
        maxval = max(img(:));
        issigned = true;
        warning('Real value arrays will be rescaled');
      case {'uint8', 'uint16', 'uint32', 'uint64'}
        minval = intmin(type);
        maxval = intmax(type);
        issigned = false;
  end

  max16 = double(intmax('uint16'));
  minval = double(minval);
  maxval = double(maxval);

  offset = -minval;
  scaling = max16 / (maxval - minval);

  infos = struct('issigned', issigned, ...
                 'offset', offset, ...
                 'scaling', scaling);

  return;
end
