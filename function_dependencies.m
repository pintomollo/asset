function [bg, alones] = functon_dependencies(path)
  
  if (nargin == 0 | isempty(path))
    path = '.';
  end

  funcs = dir([path filesep '*.m']);
  nfuncs = length(funcs);
  list = cell(nfuncs, 1);
  dependencies = cell(nfuncs, 1);
  all_names = {};

  for i = 1:nfuncs
    name = funcs(i).name;
    list{i, 1} = name;
    depends = depfun(name, '-toponly', '-quiet');
    depends = cellfun(@(c){regexprep(c,'/Appl.+','')}, depends);
    depends = cellfun(@(c){regexprep(c,'/.+/','')}, depends);

    for j=length(depends):-1:1
      if (isempty(depends{j}))
        depends(j) = [];
      elseif (strncmp(name, depends{j}, length(name)))
        depends(j) = [];
      end
    end

    dependencies{i, 1} = depends;
    all_names = [all_names; depends];
  end

  [all_names, ~, indexes] = unique([list; all_names]);

  rel_mat = zeros(length(all_names));
  for i=1:nfuncs
    
    assigns = ismember(all_names, dependencies{i, 1});
    rel_mat(indexes(i), :) = assigns;
  end

  bg = biograph(rel_mat, all_names);

  alones = all_names(~(any(rel_mat, 2)));

  return;
end
