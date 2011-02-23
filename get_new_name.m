function new_name = get_new_name(pattern, folder)

  if (findstr(pattern, '*') > 0)
    error(['Regular expression ''' pattern ''' is too vague (''*'') to get incremented']);
  end

  if (nargin == 1)
    folder = '.';
  end

  cd_dir = pwd;

  if (~exist(folder, 'dir'))
    mkdir(folder);
  end

  folder_dir = fullfile(cd_dir, folder);
  new_tokens = [];

  ls_dir = dir(folder_dir);
  for d = 1:length(ls_dir)
    [tmp tokens] = regexp(ls_dir(d).name, pattern, 'match', 'tokens');
    ntokens = length(tokens);

    if (ntokens > 0)
      tmp_tokens = zeros(1, ntokens);
      for i = 1:ntokens
        tmp_tokens(i) = str2num(char(tokens{i}));
      end
      new_tokens = [new_tokens tmp_tokens];
    end
  end

  if (length(new_tokens) > 0)
    tokens = sort(new_tokens);

    if (tokens(1) > 1)
      indxs = [1:tokens(end)+1];
    else
      indxs = [tokens(1):tokens(end)+1];
    end
    indxs = setxor(tokens, indxs);

    indx = num2str(indxs(1));
  else
    indx = '1';
  end

  new_name = regexprep(pattern, '\(.+\)', indx);
  new_name = regexprep(new_name, '[\\\?]', '');
  new_name = fullfile(folder_dir, new_name);

  'need to adapt for relative paths'
  new_name = relativepath(new_name);

  return;
end
