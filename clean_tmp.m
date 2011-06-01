function clean_tmp

  ls_dir = dir('*.mat'); 
  used_tmp = [];

  for i=1:length(ls_dir)
    load(ls_dir(i).name);

    if (~exist('mymovie', 'var'))
      continue;
    end

    fields = fieldnames(mymovie);
    for j=1:length(fields)
      if (~isempty(mymovie.(fields{j})) && isfield(mymovie.(fields{j}), 'fname'))
        for k=1:length(mymovie.(fields{j}))
          [tmp tokens] = regexp(mymovie.(fields{j})(k).fname,'tmpmat(\d+)\.*','match','tokens');
          if (length(tokens)~=0)
            used_tmp = [used_tmp; str2double(char(tokens{1}))];
          end
        end
      end
    end

    clearvars -except 'ls_dir' used_tmp;
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
        delete([tmp_dir ls_dir(d).name]);
      end
    end
  end

  return;
end
