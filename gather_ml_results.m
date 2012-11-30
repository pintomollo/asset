function results = gather_ml_results(fname, file_pattern, keep_evolution)

  if (nargin == 1)
    file_pattern = '.*%.4f?_evol\\.dat';
    keep_evolution = false;
  elseif (nargin == 2)
    if (islogical(file_pattern))
      keep_evolution = file_pattern;
      file_pattern = '.*%.4f?_evol\\.dat';
    else
      keep_evolution = false;
    end
  elseif (nargin == 3 & islogical(file_pattern))
    tmp = file_pattern;
    file_pattern = keep_evolution;
    keep_evolution = tmp;
  end

  if (isempty(file_pattern))
    file_pattern = '.*%.4f?_evol\\.dat';
  end

  if (~exist(fname, 'file'))
    results = [];
    return;
  end

  data = {};
  process = {};
  params = {};
  results = cell(0,2);

  fid = fopen(fname, 'rt');
  if (fid<0)
    fname = absolutepath(pwd, fname);
    fid = fopen(fname,'rt');
    if (fid<0)
      return;
    end
  end

  line = fgetl(fid);
  while (ischar(line))
    tmp = regexp(line, '^([\d\.]+) (\S+) OK$', 'tokens');
    line = fgetl(fid);
    if (~isempty(tmp) & ~isempty(tmp{1}))
      done = ismember(results(:,1), tmp{1}{2});
      
      if (~any(done))
        results(end+1, 1) = tmp{1}(2);
        results(end, 2) = {{tmp{1}{1}}};
      else
        results(done, 2) = {[results{done, 2}; tmp{1}(1)]};
      end
    end
    %{
    if (isempty(tmp{1}{3}))
      done = ismember(process, tmp{1}{1});
      if (sum(done) == 1)
        prev = ismember(results(:,1), data(done));
        if (any(prev))
          results(prev, 2) = {[results{prev, 2}; [process(done), params(done)]]};
        else
          results(end+1, 1) = data(done);
          results(end, 2) = {[results{end, 2}; [process(done), params(done)]]};
        end
      else
        warning(['Process ' tmp{1}{1} ' seems to terminate before it started ! Ignoring it']);
      end
      data = data(~done);
      process = process(~done);
      params = params(~done);
    else
      process = [process; tmp{1}{1}];
      data = [data; tmp{1}{2}];
      params = [params; tmp{1}{3}];
    end
    %}
  end
  fclose(fid);

  [path, pattern, ext] = fileparts(file_pattern);
  file_pattern = [pattern ext];
  if (isempty(path))
    path = pwd;
  end

  % Looking for the files
  for i=1:size(results, 1)
    data = results{i,2};
    prev_file = '';
    prev_indx = 0;

    for j=1:size(data, 1)
      num_proc = str2double(data{j, 1});

      pattern = sprintf(file_pattern, num_proc);
      files = regexpdir(path, pattern, true);
      first_res = numel(files);

      tmp_res = first_res;
      precision = 4;
      while(tmp_res ~= 1 & precision > 0)
        precision = precision - 1;
        tmp_pattern = sprintf(regexprep(file_pattern, '%\.\d', ['%.' num2str(precision)]), num_proc);
        files = regexpdir(path, tmp_pattern, true);
        tmp_res = numel(files);
      end

      if (tmp_res == 1)
        if (isempty(prev_file) | ~strcmp(prev_file, files))
          data{j, 1} = parse_ml_results(files{1}, Inf, keep_evolution, 'none');
          prev_file = files;
          prev_indx = j;
        else
          data{j, 1} = prev_indx;
        end
      elseif (first_res > 1)
        warning(['Too many files associated with pattern ' pattern]);
      else
        warning(['Cannot identify the file associated with pattern ' pattern]);
      end
    end

    %data{j+1, 1} = prev_file;
    results{i, 2} = data;
  end

  return;
end
