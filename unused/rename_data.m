function [finalname] = rename_data(fname,newname)

  new_tokens = [];

  ls_dir = dir;
  for d = 1:length(ls_dir)
    [tmp tokens] = regexp(ls_dir(d).name,[newname '(\d+)\.tmp','match','tokens');
    if(length(tokens)~=0)
      tmp_tokens = [];
      for i=1:length(tokens)
        tmp_tokens(i) = str2double(char(tokens{i}));
      end
      new_tokens = [new_tokens sort(tmp_tokens)];
    end
  end

  if(length(new_tokens)==0)
    indx = 1;
  else
    tokens = sort(new_tokens);
    indxs = [1:tokens(end)+1];
    indxs = setxor(tokens, indxs);

    indx = indxs(1);
  end
  finalname = [newname num2str(indx) '.tmp'];

  movefile(fname, finalname); 
end
