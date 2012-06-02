function [header, data, varying] = parse_mcmc_results(fname)

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
  frewind(fid);

  add_num = 15;
  add_pattern = ' %f';
  pattern = ['%f (' repmat('%f, ', 1, ncomma) '%f) :' repmat(' %f', 1, nfields)];

  data = textscan(fid, pattern);

  fclose(fid);

  data = cat(2, data{:});

  if (nargout == 1)
    header = data;
  elseif (nargout == 2)
    varying = ~all(bsxfun(@eq, data, data(1, :)), 1);
    header = data;
    data = varying;
  elseif (nargout == 3)
    header = data(:, 1:ncomma + 2);
    data = data(:, ncomma+3:end);
    varying = ~all(bsxfun(@eq, data, data(1, :)), 1);
  end

  return;
end
