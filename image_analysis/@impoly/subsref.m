## Code adapted from the one by Juan Pablo Carbajal <carbajal@ifi.uzh.ch> (2012)
## to enable the class to convert the calls "class.function()" to "function(class)"

## -*- texinfo -*-
## @deftypefn {Function File} {} function_name ()
## @end deftypefn

function varargout = subsref (obj, idx)

  persistent __method__ method4field typeNotImplemented
  if isempty(__method__)

    __method__ = struct();

    # List all the subfunctions to be converted
    __method__.getPosition = @(o,varargin) getPosition (o, varargin{:});
    __method__.setPosition = @(o,varargin) setPosition (o, varargin{:});
    __method__.getClosed = @(o,varargin) getClosed (o, varargin{:});
    __method__.setClosed = @(o,varargin) setClosed (o, varargin{:});
    __method__.delete = @(o,varargin) delete (o, varargin{:});
    __method__.display = @(o,varargin) display (o, varargin{:});

    # Error strings
    method4field = "Class #s has no field #s. Use #s() for the method.";
    typeNotImplemented = "#s no implemented for class #s.";

  end

  # Make sure the object is of the proper type
  if ( !strcmp (class (obj), 'impoly') )
    error ("Object must be of the impoly class but '#s' was used", class (obj) );
  elseif ( idx(1).type != '.' )
    error ("Invalid index for class #s", class (obj) );
  endif

  # Retrive the function handle to the proper subfunction
  method = idx(1).subs;
  if ~isfield(__method__, method)
    error('Unknown method #s.',method);
  else
    fhandle = __method__.(method);
  end

  # methods have two arguments, properties only one.
  if numel (idx) == 1 # can't access properties, only methods

    error (method4field, class (obj), method, method);

  end

  # Defines a function call
  if strcmp (idx(2).type, '()')

    args = idx(2).subs;

    # Only getPosition can return a value
    if (method(1) == 'g')
      if isempty(args)
        out = fhandle (obj);
      else
        out = fhandle (obj, args{:});
      end
      varargout{1} = out;
    else

      if isempty(args)
        fhandle (obj);
      else
        fhandle (obj, args{:});
      end
    end

  else

    error (typeNotImplemented,[method idx(2).type], class (obj));

  end

endfunction
