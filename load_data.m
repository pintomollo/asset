function [data] = load_data(fname, indexes)
%LOAD_DATA Loading from a binary temporary file
%
%    LOAD_DATA(FNAME, INDEXES) loades the frames specified in
%    INDEXES from the file FNAME. It returns them concatenated
%    as an array.

  [ssize, frames, frame_size, h] = size_data(fname);
  isflat=false;
  if(sum(ssize~=1)==1)
    isflat=true;
  end
  datasize = prod(ssize);

  indexes = indexes(indexes<=frames);
  indexes = diff([0 indexes]) - 1;

  if(isflat)
    data = zeros([datasize length(indexes)]);
  else
    data = zeros([ssize(1:2) length(indexes)]);
  end

  for i=1:length(indexes)
    fseek(h,frame_size*indexes(i),0);
    tmpimg = fread(h,datasize,'double');
    if(isflat)
      data(:,i) = reshape(tmpimg,ssize);
    else
      data(:,:,i) = reshape(tmpimg,ssize);
    end
  end
  fclose(h);

  if(isempty(data))
    data = [];
  end

  return;
end
