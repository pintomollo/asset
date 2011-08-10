function print_all(variable, spacer)

  if (nargin == 1)
    spacer = '';
  end

  if (isempty(variable))
    fprintf([spacer '[]\n']);
  else
    switch get_type(variable)
      case 'cell'
        print_cell(variable, spacer);
      case {'struct', 'errmsg'}
        print_structure(variable, spacer)
      case 'num'
        fprintf([spacer '%f\n'], variable);
      case 'bool'
        fprintf([spacer '%d\n'], variable);
      case 'char'
        fprintf([spacer '%s\n'], variable);
      otherwise
        fprintf([spacer '%s\n'], class(variable));
    end
  end
  
  return;
end

function print_structure(struct, spacer)

  spacer_unit = '    ';

  fields = fieldnames(struct);
  for i=1:numel(struct)
    if (numel(struct) > 1)
      fprintf([spacer '----%d----\n'], i);
    end
    for j=1:length(fields)
      name = fields{j};
      values = struct(i).(name);

      fprintf([spacer '%s:\n'], name);
      print_all(values, [spacer spacer_unit]);
    end
  end

  return;
end

function print_cell(variable, spacer)

  spacer_unit = '    ';

  if (numel(variable) == 1)
    print_all(variable{1}, spacer);
  else
    for i=1:numel(variable)
      print_all(variable{i}, [spacer spacer_unit]);
    end
  end

  return;
end
