function opts = load_parameters(opts, fname)
% LOAD_PARAMETERS loads parameters from a configuration file into the options structure.
%
%   OPTS = LOAD_PARAMETERS(OPTS, FNAME) loads the parameters listed in FNAME into OPTS.
%   For an example of the syntax of configuration files, see Config/example.txt.
%
%   OPTS = LOAD_PARAMETERS(OPTS) loads the parameters from OPTS.CONFIG_FILE.
%
%   OPTS = LOAD_PARAMETERS(FNAME) uses the standard options structure provided by get_struct('ASSET').
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 10.12.2010

  % In case there is only one argument, it might either be opts or fname
  if (nargin == 1)

    % If it contains this field, then it's opts
    if (isfield(opts, 'config_file'))
      fname = opts.config_file;

    % If it's a string, then it's fname
    elseif (ischar(opts))
      fname = opts;
      opts = get_struct('options',1);

    % Otherwise, we just don't do anything
    else
      return;
    end
  end

  % If the name is not provided correctly, we cannot do anything
  if (isempty(fname))
    return;
  end

  % If the file does not exists, we have a few other options
  if (~exist(fname, 'file'))

    % Maybe the extension was forgotten
    if (exist([fname '.txt'], 'file'))
      fname = [fname '.txt'];

    % Maybe it's located in the configuration folder
    elseif (exist(['Config' filesep fname], 'file'))
      fname = ['Config' filesep fname];

    % Or maybe even both previous cases
    elseif (exist(['Config' filesep fname '.txt'], 'file'))
      fname = ['Config' filesep fname '.txt'];

    % Maybe it's located in the configuration sub-folder
    elseif (exist(['celegans-analysis' filesep 'Config' filesep fname], 'file'))
      fname = ['celegans-analysis' filesep 'Config' filesep fname];

    % Or maybe even both previous cases
    elseif (exist(['celegans-analysis' filesep 'Config' filesep fname '.txt'], 'file'))
      fname = ['celegans-analysis' filesep 'Config' filesep fname '.txt'];

    % Otherwise we ran out of options
    else
      warning(['Configuration file ''' fname ''' could not be found.'])

      return;
    end
  end

  % We open it in text mode 
  fid = fopen(fname,'rt');

  % If there is an error, maybe we don't have the absolute path
  if (fid<0)
    fname = fullfile(pwd, fname);

    % And if it still does not work, then we give up
    fid = fopen(fname,'rt');
    if (fid<0)
      return;
    end
  end

  % If opts itself is some text, than we get the corresponding structure
  if (ischar(opts))
    opts = get_struct(opts);
  end

  % We can have prefixes to access subfields
  prefix = '';

  % We loop throughout the file, line by line
  line = fgetl(fid);
  while ischar(line)

    % We remove unsignificant white spaces
    line = strtrim(line);

    % We ignore empty lines
    if (length(line) > 0)

      % Prefixes start by a '#'
      if (line(1) == '#')

        % Need to avoid an error if the prefix is null
        if (length(line) > 1)

          % Extract the prefix
          prefix = line(2:end);

          % And we need a trailing dot to correctly access subfields
          if (prefix(end) ~= '.')
            prefix = [prefix '.'];
          end

        % Null prefix to reaccess fields at the root of opts
        else
          prefix = '';
        end

      % We avoid also comments which starts by '%'
      elseif (line(1) ~= '%')

        % We extract the field name (only chars) and the corresponding value
        tokens = regexp(line,'^(\w+)\s+(.+)$','tokens');

        % If we found the two elements, we can assign the value to the field
        if (~isempty(tokens) & length(tokens{1}) == 2)

          try
            % We use the eval function to interpret the values as in MATLAB 
            eval(['opts.' prefix '(tokens{1}{1}) = ' tokens{1}{2} ';']);
          catch
            break;
          end
        end
      end
    end

    % Process the next line
    line = fgetl(fid);
  end
  fclose(fid);

  opts = set_pixel_size(opts);

  return;
end
