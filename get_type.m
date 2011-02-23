function type = get_type(var)

  if (nargin < 1 || isempty(var))
    type = 'none';
  elseif (ischar(var))
    type = 'char';
  elseif (isnumeric(var))
    type = 'num';
  elseif (islogical(var))
    type = 'bool';
  elseif (iscell(var))
    type = 'cell';
  elseif (isstruct(var))
    type = 'struct';
  elseif (isjava(var))
    type = 'java';
  else
    type = 'unknown';
  end

  return;
end
