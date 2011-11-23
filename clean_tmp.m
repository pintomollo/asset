function clean_tmp

  ls_dir = dir('*.mat'); 
  used_tmp = [];

  for i=1:length(ls_dir)
    load(ls_dir(i).name);

    if (~exist('mymovie', 'var'))
      continue;
    end

    is_empty = recursive_empty(mymovie, true);

    %fields = fieldnames(mymovie);
    %for j=1:length(fields)
    %  if (~isempty(mymovie.(fields{j})))
    %    if (isfield(mymovie.(fields{j}), 'file'))
    %      for k=1:length(mymovie.(fields{j}))
    %        if (exist(mymovie.(fields{j})(k).file) == 0)
    %          is_empty = true;
    %        end
    %      end
    %    end
    %  end
    %end
    if (is_empty)
      disp(['Deleting ' ls_dir(i).name]);
      %delete(ls_dir(i).name);
    else
      tmp_tmp = recursive_tmp(mymovie);
      used_tmp = [used_tmp; tmp_tmp];

      %for j=1:length(fields)
      %  if (~isempty(mymovie.(fields{j})))
      %    if (isfield(mymovie.(fields{j}), 'fname'))
      %      for k=1:length(mymovie.(fields{j}))
      %        [tmp tokens] = regexp(mymovie.(fields{j})(k).fname,'tmpmat(\d+)\.*','match','tokens');
      %        if (length(tokens)~=0)
      %          used_tmp = [used_tmp; str2double(char(tokens{1}))];
      %        end
      %      end
      %    end
      %  end
      %end
    end

    clearvars -except 'ls_dir' 'used_tmp';
  end
  used_tmp = unique(used_tmp);

  cd_dir = [pwd filesep];
  if (exist('TmpData', 'dir'))
    tmp_dir = [cd_dir 'TmpData' filesep];
  else 
    tmp_dir = cd_dir;
  end

  ls_dir = dir([tmp_dir 'tmpmat*']);
  for d = 1:length(ls_dir)
    [tmp tokens] = regexp(ls_dir(d).name,'tmpmat(\d+)\.*','match','tokens');
    if(length(tokens)~=0)
      tmp_token = str2double(char(tokens{1}));
      if (~any(used_tmp == tmp_token))
        disp(['Deleting ' ls_dir(d).name]);
        %delete([tmp_dir ls_dir(d).name]);
      end
    end
  end

  return;
end

function is_empty = recursive_empty(mystruct, is_empty)

  fields = fieldnames(mystruct);
  for j=1:length(fields)
    if (~isempty(mystruct.(fields{j})))
      if (isfield(mystruct.(fields{j}), 'file'))
        for k=1:length(mystruct.(fields{j}))
          if (exist(mystruct.(fields{j})(k).file) ~= 0)
            is_empty = false;

            return;
          end
        end
      elseif (isstruct(mystruct.(fields{j})))
        for k=1:length(mystruct.(fields{j}))
          is_empty = recursive_empty(mystruct.(fields{j})(k), is_empty);

          if (~is_empty)
            return;
          end
        end
      end
    end
  end

  return;
end

function tmps = recursive_tmp(mystruct)

  tmps = [];
  used_tmp = [];

  fields = fieldnames(mystruct);
  for j=1:length(fields)
    if (~isempty(mystruct.(fields{j})))
      if (isfield(mystruct.(fields{j}), 'fname'))
        for k=1:length(mystruct.(fields{j}))
          [tmp tokens] = regexp(mystruct.(fields{j})(k).fname,'tmpmat(\d+)\.*','match','tokens');
          if (length(tokens)~=0)
            used_tmp = [used_tmp; str2double(char(tokens{1}))];
          end
        end
      elseif (isstruct(mystruct.(fields{j})))
        for k=1:length(mystruct.(fields{j}))
          tmp_tmp = recursive_tmp(mystruct.(fields{j})(k));

          tmps = [tmps; tmp_tmp];
        end
      end
    end
  end

  tmps = [tmps; used_tmp];

  return;
end
