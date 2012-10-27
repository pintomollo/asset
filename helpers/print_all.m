function print_all(varargin)

  fid = 1;
  spacer = '';
  variable = [];
  prefix = '';

  for i=1:length(varargin)
    if (isnumeric(varargin{i}) & ~isempty(fopen(varargin{i})))
      fid = varargin{i};
    elseif (isempty(variable))
      variable = varargin{i};
    else
      spacer = varargin{i};
    end
  end

  if (isempty(spacer))
    if (fid > 2)
      spacer = '\t';
    else
      spacer = ': ';
    end
  end

  myprint(fid, variable, spacer, prefix);

  return;
end

function myprint(fid, variable, spacer, prefix)

  orig_prefix = prefix;

  if (isempty(variable))
    fprintf(fid, [prefix spacer '[]\n']);
  else
    switch class(variable)
      case 'cell'
        for i=1:numel(variable)
          myprint(fid, variable{i}, spacer, [prefix '{' num2str(i) '}']);
        end
      case {'struct', 'MException'}
        fields = fieldnames(variable);
        for i=1:numel(variable)
          if (numel(variable) > 1)
            prefix = [orig_prefix '(' num2str(i) ')'];
          end
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
      case {'double', 'single', 'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32', 'uint64'}
        if (numel(variable) == 1)
          fprintf(fid, [prefix spacer '%f\n'], variable);
        else
          fprintf(fid, [prefix spacer '[']);
          fprintf(fid, '%f ', variable);
          fprintf(fid, ']\n');
        end
      case 'logical'
        if (numel(variable) == 1)
          fprintf(fid, [prefix spacer '%d\n'], variable);
        else
          fprintf(fid, [prefix spacer '[']);
          fprintf(fid, '%d ', variable);
          fprintf(fid, ']\n');
        end
      case 'char'
        fprintf(fid, [prefix spacer '''%s''\n'], variable);
      otherwise
        fprintf(fid, [prefix spacer '''%s''\n'], class(variable));
    end
  end
  
  return;
end
