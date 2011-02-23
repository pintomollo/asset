function [fname] = modify_data(fname, data, indexes)

  if(isempty(fname) | (isstruct(fname) & isfield(fname, 'fname') & isempty(fname.fname)))
    fname = store_data(fname, data);
    return;
  end

  [ssize, nframes, frame_size, h] = size_data(fname);

  if(length(indexes)==1 && indexes>nframes)
    fname = store_data(fname, data);
    return;
  end

  datasize = prod(ssize);
  
  indexes = indexes(indexes<=nframes);
  indexes = diff([0 indexes]) - 1;

  if(isfield(data,'cdata'))
    for i=1:length(indexes)
      fseek(h,frame_size*indexes(i),0);
      count = fwrite(h,data(i).cdata,'double');
    end
  else
    for i=1:length(indexes)
      fseek(h,frame_size*indexes(i),0);
      count = fwrite(h,squeeze(data(:,:,i)),'double');
    end
  end
  fclose(h);
    
  return;
end
