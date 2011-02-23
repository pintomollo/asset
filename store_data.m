function [fname] = store_data(fname, data)

  if (isstruct(fname) & isfield(fname, 'fname'))
    h = fname.fname;
  else
    h = fname;
  end

  if (prod(size(data))==0)
    return
  end

  if(~isinteger(h))
    if(isempty(h))
      
      cd_dir = cd;
      if (exist('TmpData', 'dir'))
        tmp_dir = [cd_dir '/TmpData/'];
      elseif (exist('../TmpData','dir'))
        cd('../TmpData/');
        tmp_dir = cd;
        cd(cd_dir);
      else 
        tmp_dir = cd_dir;
      end

      new_tokens = [];

      ls_dir = dir(tmp_dir);
      for d = 1:length(ls_dir)
        [tmp tokens] = regexp(ls_dir(d).name,'tmpmat(\d+)\.tmp','match','tokens');
        if(length(tokens)~=0)
          tmp_tokens = [];
          for i=1:length(tokens)
            tmp_tokens(i) = str2num(char(tokens{i}));
          end
          new_tokens = [new_tokens sort(tmp_tokens)];
        end
      end

      if (length(new_tokens) > 0)
        tokens = sort(new_tokens);
        indxs = [1:tokens(end)+1];
        indxs = setxor(tokens, indxs);

        indx = indxs(1);
      else
        indx = 1;
      end

      new_name = fullfile(tmp_dir, ['tmpmat' num2str(indx) '.tmp']);
      new_name = relativepath(new_name);

      h = fopen(new_name,'a');

      if(isfield(data,'cdata'))
        ssize = size(data(1).cdata);
      elseif(isfield(data,'data'))
        ssize = size(data(1).data);
      else
        ssize = size(data);
      end
      if(length(ssize)<3)
        ssize = [ssize ones(1,3-length(ssize))];
      end
      count = fwrite(h,ssize(1:3),'double');

    elseif(ischar(h))
      new_name = h;
      h = fopen(new_name,'a');
    end
  else
    new_name = h;
    frewind(h);
  end

  if(isfield(data,'cdata'))
    for i=1:length(data)
      count = fwrite(h,data(i).cdata,'double');
    end
  elseif(isfield(data,'data'))
    for i=1:length(data)
      count = fwrite(h,data(i).data,'double');
    end
  else
    count = fwrite(h,data,'double');
  end

  if (isstruct(fname) & isfield(fname, 'fname'))
    fname.fname = new_name;
  else
    fname = new_name;
  end
  
  fclose(h);
end
