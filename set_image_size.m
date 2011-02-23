function params = set_image_size(params, imgsize)

  fields = fieldnames(params);

  for i=1:length(fields)
    if (ischar(params.(fields{i})))
      commands = params.(fields{i});
      indx = [0 strfind(commands, ';')];

      for j=1:length(indx)-1
        if (j == length(indx) - 1)
          eval(['params.(fields{i}) = ' commands(indx(j)+1:indx(j+1))]);
        else
          eval(commands(indx(j)+1:indx(j+1)));
        end
      end
    elseif (isstruct(params.(fields{i})))
      params.(fields{i}) = set_image_size(params.(fields{i}), imgsize);
    end
  end

  return;
end
