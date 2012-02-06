function type = get_type(var)
% GET_TYPE returns a string defining the type of the input variable. This function
% simply combines the different build-in function (ischar, isnumeric, ...) into a
% string which is very handy in switch statements.
%
%   TYPE = GET_TYPE(VARIABLE) returns the type of VARIABLE as a string. Possible
%   outputs are the following:
%     'none'    No type, either no input or an empty variable
%     'char'    String input      help ischar
%     'num'     Numerical input   help isnumeric
%     'bool'    Logical input     help islogical
%     'cell'    Cell input        help iscell
%     'struct'  Structure input   help isstruct
%     'java'    Java object       help isjava
%     'handle'  Function handle   help function_handle
%     'errmsg'  Matlab error      help MException
%     'unkown'  Unkown type, none of the above
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 06.02.2012

  % Simply check the different possible types, one-by-one, using the build-in
  % Matlab functions.
  if (nargin < 1 || isempty(var))
    type = 'none';
  elseif (isnumeric(var))
    type = 'num';
  else
    type = class(var);
    switch type
      case {'char', 'cell', 'struct'}
      case 'logical'
        type = 'bool';
      case 'function_handle'
        type = 'handle';
      otherwise
        if (isjava(var))
          type = 'java';
        elseif (strncmp(type, 'MException', 10))
          type = 'errmsg';
        else
          type = 'unknown';
        end
    end
  end

  return;
end
