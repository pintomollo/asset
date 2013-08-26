function sort_ml_results(varargin)

  [fnames, groups, filters, folder, type] = parse_input(varargin{:});

  [path, name, ext] = fileparts(fnames);
  new_folder = fullfile(path, folder);

  if (mkdir(new_folder))
    datas = group_ml_results(fnames, groups, filters);
    selected_files = cell(size(datas,1), 1);

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

          selected_files{i} = vals{best_indx, 1};
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

function [fnames, groups, filters, folder, type] = parse_input(varargin)

  fnames = 'adr-kymo-*_evol.dat';
  groups = {'type'; 'parameter_set'};
  filters = {};
  folder = '';
  type = 'best';

  for i=1:nargin
    var_type = class(varargin{i});
    switch var_type
      case 'char'
        name = varargin{i};
        if (any(name == '*'))
          fnames = name;
        elseif (isempty(folder))
          folder = name;
        else
          type = name;
        end
      case 'cell'
        if (size(varargin{i},2) == 1)
          groups = varargin{i};
        else
          filters = varargin{i};
        end
    end
  end

  if (isempty(folder))
    folder = 'Sorted';
  end

  return;
end
