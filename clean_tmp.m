function clean_tmp(do_it_really)

  if (nargin == 0)
    do_it_really = true;
  end

  ls_dir = dir('*.mat'); 
  used_tmp = [];

  for i=1:length(ls_dir)
    disp(['Checking in ' ls_dir(i).name '...']);

    load(ls_dir(i).name);

    if (~exist('mymovie', 'var'))
      continue;
    end
    %is_empty = recursive_empty(mymovie);

    %if (is_empty)
    %  disp(['Deleting ' ls_dir(i).name]);
%
%      if (do_it_really)
%        delete(ls_dir(i).name);
%      end
%    else
      tmp_tmp = recursive_tmp(mymovie);
      used_tmp = [used_tmp; tmp_tmp];
%    end

    clearvars -except 'ls_dir' 'used_tmp' 'do_it_really';
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

        if (do_it_really)
          delete([tmp_dir ls_dir(d).name]);
        end
      end
    end
  end

  return;
end

function is_empty = recursive_empty(mystruct)

  fields = fieldnames(mystruct);

  indx = ismember(fields, 'file');
  is_empty = true;

  if (any(indx))
    if (exist(mystruct.file) ~= 0)
      is_empty = false;

      return;
    end

    fields = fields(~indx);
  end

  for j=1:length(fields)
    if (isstruct(mystruct.(fields{j})))
      for k=1:length(mystruct.(fields{j}))
        if (k==1)
          sub_fields = fieldnames(mystruct.(fields{j})(k));
          has_struct = any(ismember(sub_fields, 'file'));
          if (~has_struct)
            for s = 1:length(sub_fields)
              if (isstruct(mystruct.(fields{j})(k).(sub_fields{s})))
                has_struct = true;
                break;
              end
            end
          end
        end

        if (has_struct)
          is_empty = recursive_empty(mystruct.(fields{j})(k));

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

  indx = ismember(fields, 'fname');

  if (any(indx))
    [tmp tokens] = regexp(mystruct.fname,'tmpmat(\d+)\.*','match','tokens');
    if (length(tokens)~=0)
      used_tmp = [used_tmp; str2double(char(tokens{1}))];

      if (~exist(mystruct.fname, 'file'))
        diaply(['Missing ' mystruct.fname]);
      end
    end

    fields = fields(~indx);
  end

  for j=1:length(fields)
    if (isstruct(mystruct.(fields{j})))
      for k=1:length(mystruct.(fields{j}))
        if (k==1)
          sub_fields = fieldnames(mystruct.(fields{j})(k));
          has_struct = any(ismember(sub_fields, 'file'));
          if (~has_struct)
            for s = 1:length(sub_fields)
              if (isstruct(mystruct.(fields{j})(k).(sub_fields{s})))
                has_struct = true;
                break;
              end
            end
          end
        end

        if (has_struct)
          tmp_tmp = recursive_tmp(mystruct.(fields{j})(k));

          tmps = [tmps; tmp_tmp];
        end
      end
    end
  end

  tmps = [tmps; used_tmp];

  return;
end
