function params = set_pixel_size(params, pixel_size)
% SET_PIXEL_SIZE computes the actual size of the pixel in the image using the CCD
% camera pixel size, the magnification and the binning. In addition, it computes the
% other values in the parameter structure which values depend on "pixel_size"
% (see get_struct.m).
% 
%   OPTS = SET_PIXEL_SIZE(OPTS) computes pixel_size using the fields 'ccd_pixel_size',
%   'magnification' and 'binning' of the parameter structure OPTS (get_struct('options')).
%   Using this value, it then computes the actual value of the dynamic fields of OPTS.
%
%   OPTS = SET_PIXEL_SIZE(OPTS, PIXEL_SIZE) uses the provided value for pixel_size
%   to set the value of the dynamic fields.
%
%   OPTS = SET_PIXEL_SIZE uses the default value for OPTS too.
%
% NOTES ON USING DYNAMIC FIELD VALUES
%   In case the value of the parameter you want to set depends on the size of the pixel
%   (in um), you can substitute its value by the expression that computes it. In this
%   expression, you can then use the variable 'pixel_size' which will contain the
%   appropriate value. Simply put the whole expression in a string terminated by the
%   usual ';' character. This function will then evaluate this expression using the
%   value of 'pixel_size' and store the result into the same field.
%
%   Multiple lines can be stored in a single field separated by ';'. This technique
%   allows to define your own temporary variables. In this case you will need to
%   explicit the assignment as in any normal expression using '='.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 10.12.2010

  % Check if we need to get a new option structure
  if (nargin == 0)
    params = get_struct('options');
    params = set_pixel_size(params);

    return;

  % Otherwise, we have only the option structure provided
  elseif (nargin == 1 | isempty(pixel_size))

    % We need at least these two fields to exist
    if (isfield(params, 'ccd_pixel_size') & isfield(params, 'magnification'))

      % Maybe there is even binning
      if (isfield(params, 'binning'))
        pixel_size = params.binning * params.ccd_pixel_size / params.magnification;
      else
        pixel_size = params.ccd_pixel_size / params.magnification;
      end

      % Store the pixel size
      params.pixel_size = pixel_size;
    else
      pixel_size = [];
    end
  end

  % If we have no data for the pixel_size, stop here
  if (isempty(params))
    return;
  end

  % Otherwise, we loop recursively on the parameter structure to try to set the
  % dynamic fields
  fields = fieldnames(params);
  for n=1:length(params)
    for i=1:length(fields)

      % If it's a string, maybe we can set something here
      if (ischar(params(n).(fields{i})) & ~isempty(pixel_size))

        % Get the comamnds and split them per line
        commands = params(n).(fields{i});
        indx = [0 strfind(commands, ';')];
        ncomma = length(indx) - 1;

        % Loop over each command and execute it
        for j=1:ncomma

          % For the last one, we need to assign it back the structure field
          if (j == ncomma)
            eval(['params(n).(fields{i}) = ' commands(indx(j)+1:indx(j+1))]);
          else
            eval(commands(indx(j)+1:indx(j+1)));
          end
        end

      % If it's a structure field, call set_pixel_size recursively
      elseif (isstruct(params(n).(fields{i})))
        params(n).(fields{i}) = set_pixel_size(params(n).(fields{i}), pixel_size);
      end
    end
  end

  return;
end
