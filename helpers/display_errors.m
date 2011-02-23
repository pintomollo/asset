function [herr, stds, groups, names] = display_errors(trackings, indx, colors, subindx)
  
  if (nargin < 2)
    indx = zeros(0,3);
    colors = zeros(0,3);
    subindx = [];
  elseif (nargin < 3)
    colors = zeros(0,3);
    subindx =[];
  elseif (nargin < 4)
    subindx = [];
  end

  if (isempty(indx))
    indx = zeros(1,3);
  end

  [errs, stds, groups, names] = extract_errors(trackings, indx);

  if (isempty(colors))
    colors = hsv(6* ceil(size(errs,2) / 6));
  end

  if (~isempty(subindx))
    errs = errs(:, subindx);
    stds = stds(:, subindx);
    names = names(subindx);
  end

  if (nargout > 1)
    herr = errs;
  elseif (~isempty(errs))
    figure;
    herr = barweb(errs * 100, stds * 100, [], groups, 'Average distance to the mean tracking', 'Mean tracking', 'Average distance (% of embryo''s radius)', colors, [], names);
  end
end

function [errs, stds, groups, names] = extract_errors(trackings, indx)

  errs = [];
  stds = [];
  groups = {};
  names = {};

  if (isstruct(trackings))
    nchild = 0;
    nfiles = 0;

    if (isfield(trackings, 'child'))
      nchild = length(trackings.child);
    end
    
    for  i=1:nchild 
      [child_errs, child_stds, child_groups, child_names] = extract_errors(trackings.child(i), indx);

      if (~isempty(child_errs) && ~all(child_errs(:) == 0))
        errs = [errs child_errs];
        stds = [stds child_stds];
        names = [names child_names];

        if (size(errs, 1) > length(groups))
          groups = [groups child_groups];
        end
      end
    end

    if (~isempty(errs))
      errs = [errs mymean(errs,2)];
      stds = [stds mymean(stds, 2)];
      names = [names {'Intra'}];
    end

    if (isfield(trackings, 'errors') && ~isempty(trackings.errors))
      res = trackings.errors;

      eval_indx = cell(size(indx));
      eval_indx(indx == 0) = {':'};
      eval_indx(indx ~= 0) = regexp(num2str(indx(indx ~=0)), '\w+', 'match');

      eval(['res = res(' eval_indx{1} ',' eval_indx{2} ',' eval_indx{3} ',:);']);
      res = mean(res, 4);
      res = reshape(res, size(res, 1), []);

      [mres, sres] = mymean(res, 2);

      errs = [errs mres];
      stds = [stds sres];

      if (isfield(trackings, 'name') && ~isempty(trackings.name))
        names = [names trackings.name];
      else
        names = [names {'Inter'}];
      end
    end

    if (isfield(trackings, 'files') && length(groups) < size(errs, 1))
      nfiles = length(trackings.files);

      for i=1:nfiles
        if (isfield(trackings.files(i), 'groups'))
          groups = [groups trackings.files(i).groups];
        end

        if (length(groups) > size(errs, 1) && indx(1) ~= 0)
          groups = groups(indx(1));
        end

        if (length(groups) >= size(errs, 1))
          break;
        end
      end
    end
  end

  return;
end
