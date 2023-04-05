function save_parameters(opts, fname)
% SAVE_PARAMETERS saves the content of a parameter structure (or any other structure)
% in a text file, which can then be reloaded using LOAD_PARAMETERS.
%
%   SAVE_PARAMETERS(OPTS) saves OPTS into a text file provided by the user through
%   a GUI.
%
%   SAVE_PARAMETERS(OPTS, FNAME) saves OPTS into FNAME.
%
%   SAVE_PARAMETERS(OPTS, FID) saves OPTS using the file identifier FID.
%
%   Note that SAVE_PARAMETERS alway overwrites the content of an existing file.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 14.05.14

  % First check if we need to ask for the target file
  if (nargin == 0)
    return;
  elseif (nargin < 2)
    % If so, check if we can save them in the Config folder
    if (exist('Config', 'dir'))
      conf_dir = which('Config.');
    elseif (exist(['cell-tracking' filesep 'Config'], 'dir'))
      conf_dir = ['cell-tracking' filesep 'Config'];
    else
      conf_dir = pwd;
    end

    % Ask the user for a filename
    [fname, pathname] = uiputfile({'*.txt','All text files'; '*.*','All files' }, ...
                                   'Save parameters', [conf_dir filesep]);

    % This means the user has canceled
    if (all(fname == 0))
      return;
    end

    % Get the full name
    fname = fullfile(pathname, fname);
  end

  % Check if we got a file ID already
  if (isnumeric(fname) && ~isempty(fopen(fname)))
    fid = fname;
  else
    % We open it in text mode 
    fid = fopen(fname,'w+t');

    % If there is an error, maybe we don't have the absolute path
    if (fid<0)
      fname = fullfile(pwd, fname);

      % And if it still does not work, then we skip this file
      fid = fopen(fname,'w+t');
      if (fid<0)
        return;
      end
    end
  end

  % Define the spacer and the inital prefix
  spacer = '\t';
  prefix = '';

  % Some fancy header for the file. The % indicates a comment.
  fprintf(fid, '%%Parameters saved on %s\n\n', datestr(now));

  % Try to perform the actual writing
  try
    myprint(fid, opts, spacer, prefix);
  catch ME
    warning('ASSET:save_parameters', ['An error occured when saving the parameters:\n' ME.message])
  end

  % If we did not get a file identifier, we need to close the file
  if (any(fid ~= fname))
    fclose(fid);
  end

  return;
end

function myprint(fid, variable, spacer, prefix)
% This function performs the actual writing of the structure into the file. This
% second function is necessary as structures are written recursively.

  % Keep track of the prefix
  orig_prefix = prefix;

  % Empty variables are a bit different that the others
  if (isempty(variable))
    fprintf(fid, [prefix spacer '[]\n']);
  else

    % Each type of variable deserves a different printing
    switch class(variable)

      % To simplify things, we print cell arrays in a linear fashion,
      % we hard-code the corresponding index into the prefix and we use
      % myprint recursively to print the content of each cell.
      case 'cell'
        for i=1:numel(variable)
          myprint(fid, variable{i}, spacer, [prefix '{' num2str(i) '}']);
        end

      % For the structures, we run through the fields and call recursively
      % this function to print the content of each field.
      case {'struct', 'MException'}

        % In addition to the different fields, structures can be arrays.
        fields = fieldnames(variable);
        for i=1:numel(variable)

          % Modify the prefix in case of an array
          if (numel(variable) > 1)
            prefix = [orig_prefix '(' num2str(i) ')'];
          end

          % Loop through the fields and call myprint again, adapting the prefix
          for j=1:length(fields)
            name = fields{j};
            values = variable(i).(name);

            if (isempty(prefix))
              myprint(fid, values, spacer, name);
            else
              myprint(fid, values, spacer, [prefix '.' name]);
            end
          end
        end
        % Only fancying here !
        fprintf(fid, '\n');

      % All the numbers are printed with scientific notation for maximal precision
      case {'double', 'single', 'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32', 'uint64'}
        if (numel(variable) == 1)
          fprintf(fid, [prefix spacer '%e\n'], variable);
        else
          fprintf(fid, [prefix spacer '[']);
          fprintf(fid, '%e ', variable);
          fprintf(fid, ']\n');
        end

      % Print logicals as binary values
      case 'logical'
        if (numel(variable) == 1)
          fprintf(fid, [prefix spacer '%d\n'], variable);
        else
          fprintf(fid, [prefix spacer '[']);
          fprintf(fid, '%d ', variable);
          fprintf(fid, ']\n');
        end

      % Print the strings along with the required apostrophes
      case 'char'
        fprintf(fid, [prefix spacer '''%s''\n'], variable);

      % Print the functions as handlers
      case 'function_handle'

        % For loading, we need the @
        string = func2str(variable);
        if (string(1)~='@')
          string = ['@' string];
        end

        fprintf(fid, [prefix spacer '%s\n'], string);

      % When we do not know what to do, we print the type of object we found
      otherwise
        fprintf(fid, [prefix spacer '''%s''\n'], class(variable));
    end
  end

  return;
end
