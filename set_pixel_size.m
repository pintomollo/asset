function params = set_pixel_size(params, pixel_size)
% SET_PIXEL_SIZE computes the actual size of the pixel in the image using the CCD
% camera pixel size, the magnification and the binning. In addition, it computes the
% other values in the parameter structure which values depend on "pixel_size"
% (see get_struct.m).
% 
%   OPTS = SET_PIXEL_SIZE(OPTS) computes pixel_size using the fields 'ccd_pixel_size'
%   and 'magnification' of a generic parameter structure OPTS (e.g. get_struct('ASSET')).
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
%   See "get_struct.m" for some examples of use.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 10.12.2010

  if (nargin == 0)
    params = get_struct('ASSET');
    params = set_pixel_size(params);

    return;
  elseif (nargin == 1 | isempty(pixel_size))
    if (isfield(params, 'ccd_pixel_size') & isfield(params, 'magnification'))
      if (isfield(params, 'binning'))
        pixel_size = params.binning * params.ccd_pixel_size / params.magnification;
      else
        pixel_size = params.ccd_pixel_size / params.magnification;
      end
      params.pixel_size = pixel_size;
    else
      pixel_size = [];
    end
  end

  if (isempty(params))
    return;
  end

  fields = fieldnames(params);

  for i=1:length(fields)
    if (ischar(params.(fields{i})) & ~isempty(pixel_size))
      commands = params.(fields{i});
      indx = [0 strfind(commands, ';')];
      ncomma = length(indx) - 1;

      for j=1:ncomma
        if (j == ncomma)
          eval(['params.(fields{i}) = ' commands(indx(j)+1:indx(j+1))]);
        else
          eval(commands(indx(j)+1:indx(j+1)));
        end
      end
    elseif (isstruct(params.(fields{i})))
      params.(fields{i}) = set_pixel_size(params.(fields{i}), pixel_size);
    end
  end

  return;
end
