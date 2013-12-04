function sort_ml_results(varargin)

  [fnames, groups, filters, folder, type, args] = parse_input(varargin{:});

  [path, name, ext] = fileparts(fnames);
  new_folder = fullfile(path, folder);

  if (mkdir(new_folder))
    datas = group_ml_results(fnames, groups, filters);
    selected_files = {};

    switch type
      case 'best'
        for i=1:size(datas,1)
          vals = datas{i,2};
          best_score = vals{1,2}(end).score;
          best_indx = 1;

          for j=2:size(vals,1)
            if (vals{j,2}(end).score < best_score)
              best_score = vals{j,2}(end).score;
              best_indx = j;
            end
          end

          selected_files{end+1} = vals{best_indx, 1};
        end

      case '~nparams'
        for i=1:size(datas,1)
          vals = datas{i,2};
          for j=1:size(vals,1)
            if (length(vals{j,2}(end).params) ~= args{1}(1))
              selected_files{end+1} = vals{j, 1};
            end
          end
        end

      otherwise
        warning(['Sorting type ' type ' has not yet been implemented.']);
    end
  end

  for i=1:length(selected_files)
    if (~isempty(selected_files{i}))
      copyfile(selected_files{i}, new_folder);
    end
  end

  return;
end

function [fnames, groups, filters, folder, type, args] = parse_input(varargin)

  fnames = 'adr-kymo-*_evol.dat';
  groups = {'type'; 'parameter_set';'fit_model';'fit_flow'};
  filters = {};
  folder = 'Sorted';
  args = {};
  type = '';

  for i=1:nargin
    var_type = class(varargin{i});
    switch var_type
      case 'char'
        name = varargin{i};
        if (any(name == '*'))
          fnames = name;
        elseif (isempty(type))
          type = name;
        else
          folder = name;
        end
      case 'cell'
        if (size(varargin{i},2) == 1)
          groups = varargin{i};
        else
          filters = varargin{i};
        end
      otherwise
        args{end+1} = varargin{i};
    end
  end

  if (isempty(type))
    type = 'best';
  end

  return;
end
