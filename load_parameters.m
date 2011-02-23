function opts = load_parameters(fname, opts)
  
  if (nargin == 1)
    opts = fname;

    if (isfield(opts, 'config_file'))
      fname = opts.config_file;
    elseif (ischar(fname))
      opts = get_struct('RECOS',1);
    else
      return;
    end
  end

  if (isempty(fname))
    return;
  end

  if (~exist(fname, 'file'))
    if (exist([fname '.txt'], 'file'))
      fname = [fname '.txt'];
    elseif (exist(['Config/' fname], 'file'))
      fname = ['Config/' fname];
    elseif (exist(['Config/' fname '.txt'], 'file'))
      fname = ['Config/' fname '.txt'];
    elseif (exist(['../Config/' fname], 'file'))
      fname = ['../Config/' fname];
    elseif (exist(['../Config/' fname '.txt'], 'file'))
      fname = ['../Config/' fname '.txt'];
    else
      return;
    end
  end

  fid = fopen(fname,'rt');
  if (fid<0)
    fname = [pwd fname];
    fid = fopen(fname,'rt');
    if (fid<0)
      return;
    end
  end

  prefix = '';
  line = fgetl(fid);
  while ischar(line)
    if (length(line) > 0 && line(1) == '#')
      if (length(line) > 1)
        prefix = strtrim(line(2:end));

        if (prefix(end) ~= '.')
          prefix = [prefix '.'];
        end
      else
        prefix = '';
      end
    elseif (length(line) > 0 && line(1) ~= '%')
    
      tokens = regexp(line,'^\s*(\w+)\s+(.+)\s*$','tokens');

      if (length(tokens{1}) == 2)
        field_type = identify_field_type(opts, prefix, tokens{1}{1});

        switch field_type
          case 'char'
            eval(['opts.' prefix '(tokens{1}{1}) = tokens{1}{2};']);
          case 'num'
            eval(['opts.' prefix '(tokens{1}{1}) = str2double(tokens{1}{2});']);
          case 'bool'
            eval(['opts.' prefix '(tokens{1}{1}) = logical(' tokens{1}{2} ');']);
          case 'cell'
            eval(['opts.' prefix '(tokens{1}{1}) = regexp(tokens{1}{2}, ''\W'', ''split'');']);
        end
      end
    end

    line = fgetl(fid);
  end

  return;
end

function type = identify_field_type(mystruct, prefix, field)

  indx = strfind(prefix, '.');
  if (~isempty(indx))
    substruct = prefix(1:indx(1)-1);

    if (isfield(mystruct, substruct))
      type = identify_field_type(mystruct.(substruct), prefix(indx(1)+1:end), field);
    else
      type = 'none';
    end
  else
    type = get_type(mystruct.(field));
  end

  return;
end
