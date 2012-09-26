function [header, data, varying, opts] = parse_mcmc_results(fname)

  data = [];

  fid = fopen(fname, 'rt');
  if (fid<0)
    fname = absolutepath(fname);
    fid = fopen(fname,'rt');
    if (fid<0)
      return;
    end
  end

  %{
  line = fgetl(fid);
  pattern = '';
  while ischar(line)
    if (length(line) > 0)
      if (isempty(pattern))
        indx = findstr(line, ' ');
      end
    end
    line = fgetl(fid);
    nline = nline + 1;
  end
  %}

  line = fgetl(fid);
  ncomma = length(findstr(line, ','));
  nfields = length(findstr(line, ' ')) - ncomma - 2;

  has_header = (ncomma == 0);

  nbits = 0;
  cur = ftell(fid);

  while (ncomma == 0 & ischar(line))
    line = fgetl(fid);
    ncomma = length(findstr(line, ','));

    nbits = cur;
    cur = ftell(fid);
  end
  %frewind(fid);
  fseek(fid, nbits - cur, 0);

  ncomma = length(findstr(line, ','));
  nfields = length(findstr(line, ' ')) - ncomma - 2;

  add_num = 15;
  add_pattern = ' %f';
  pattern = ['%f (' repmat('%f, ', 1, ncomma) '%f) :' repmat(' %f', 1, nfields)];

  %data = fscanf(fid, pattern)
  data = textscan(fid, pattern);

  fclose(fid);

  data = cat(2, data{:});

  if (has_header)
    opts = load_parameters(get_struct('fitting'), fname);
  else
    opts = [];
  end

  if (nargout == 1)
    header = data;
  elseif (nargout == 2)
    if (has_header)
      header = data;
      data = opts;
    else
      varying = ~all(bsxfun(@eq, data, data(1, :)), 1);
      header = data;
      data = varying;
    end
  else
    header = data(:, 1:ncomma + 2);
    data = data(:, ncomma+3:end);
    varying = ~all(bsxfun(@eq, data, data(1, :)), 1);
    if (has_header & nargout == 3)
      varying = opts;
    end
  end

  return;
end
