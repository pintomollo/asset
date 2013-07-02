function datas = group_ml_results(fnames, groups, filters, ref_file)

  datas = {};

  if (nargin == 0)
    ref_file = '';
    return;
  elseif (nargin == 1)
    groups = {};
    filters = {};
    ref_file = '';
  elseif (nargin == 2)
    if (ischar(groups))
      ref_file = groups;
      groups = {};
      filters = {};
    elseif (size(groups, 2) > 1)
      filters = groups;
      groups = {};
      ref_file = '';
    else
      filters = {};
      ref_file = '';
    end
  elseif (nargin == 3)
    tmp_ref = '';
    tmp_group = {};
    tmp_filt = {};

    if (ischar(groups))
      tmp_ref = groups;
    elseif (size(groups, 2) > 1)
      tmp_filt = groups;
    else
      tmp_group = groups;
    end
 
    if (ischar(filters))
      tmp_ref = filters;
    elseif (size(filters, 2) == 1)
      tmp_group = filters;
    else
      tmp_filt = filters;
    end

    [groups, filters, ref_file] = deal(tmp_group, tmp_filt, tmp_ref);
  end

  ngroups = size(groups, 1);
  nfilters = size(filters, 1);

  [path, pattern, ext] = fileparts(fnames);
  file_pattern = [pattern ext];
  if (isempty(path))
    path = pwd;
  end

  if (isempty(ref_file))
    files = dir(fnames);

    fields = fieldnames(files);
    files = struct2cell(files);

    files = files(ismember(fields, 'name'), :).';
  else
    files = cell(0,1);
    file_pattern = 'adr-kymo-%.4f?_evol\\.dat';

    fid = fopen(ref_file, 'rt');
    if (fid<0)
      fname = absolutepath(pwd, ref_file);
      fid = fopen(ref_file,'rt');
      if (fid<0)
        return;
      end
    end

    line = fgetl(fid);
    while (ischar(line))
      tmp = regexp(line, '^([\d\.]+) (\S+) OK$', 'tokens');
      line = fgetl(fid);
      if (~isempty(tmp) & ~isempty(tmp{1}))
        %files(end+1, 1) = tmp{1}(1);
        tmp_pattern = sprintf(regexprep(file_pattern, '%\.\d', '%.4'), str2double(tmp{1}{1}));
        tmp_file = regexpdir(path, tmp_pattern, false);
        if (~isempty(tmp_file))
          files{end+1, 1} = relativepath(tmp_file{1});
        end
      end
    end
    fclose(fid);
  end

  nfiles = length(files);
  datas = cell(nfiles, 2);
  nstored = 0;

  display('Parsing files...');

  for i=1:nfiles
    fname = [path filesep files{i}];
    opt = load_parameters(get_struct('fitting'), fname);

    is_valid = true;
    for j=1:nfilters
      is_valid = is_valid && isfield(opt, filters{j, 1}) && numel(opt.(filters{j,1})) == numel(filters{j,2}) && all(opt.(filters{j,1})(:) == filters{j,2}(:));

      if (~is_valid)
        break;
      end
    end

    for j=1:ngroups
      is_valid = is_valid && isfield(opt, groups{j, 1});

      if (~is_valid)
        break;
      end
    end

    if (is_valid)
      found = false;
      for n=1:nstored
        ref = datas{n, 1}{1};
        all_same = true;
        for j=1:ngroups
          all_same = all_same && numel(ref.(groups{j})) == numel(opt.(groups{j})) && all(ref.(groups{j})(:) == opt.(groups{j})(:));
        end

        if (all_same)
          found = n;
          break;
        end
      end

      if (~found)
        found = nstored + 1;
        datas{found, 1} = {opt};
        datas{found, 2} = {fname};
        nstored = nstored + 1;
      else
        datas{found, 2} = [datas{found, 2}; fname];
      end

      display(files{i})
    end
  end

  display('Loading files...');

  for i=1:nstored
    indx = 0;
    files = datas{i, 2};
    for j=1:length(files)
      display(files{j});
      tmp = parse_ml_results(files{j}, Inf, true, 'none');
      if (isfinite(tmp(end).score))
        indx = indx + 1;
        datas{i, 2}{indx,2} = tmp;
        datas{i, 2}{indx,1} = datas{i, 2}{j,1};
      end
    end
    datas{i,2} = datas{i,2}(1:indx,:);
  end

  datas = datas(1:nstored, :);

  return;
end
