function [data] = load_slice(fname, indexes)
%LOAD_SLICE Loading from a binary temporary file
%
%    LOAD_SLICE(FNAME, INDEXES) loades the slices specified in
%    INDEXES from the file FNAME. It returns them concatenated
%    as an array. A slice is the concatenation of the same row
%    in every frame of the stack

  [ssize, frames] = size_data(fname);
  isflat=false;
  if(sum(ssize~=1)==1)
    isflat=true;
  end

  indexes = indexes(indexes<=ssize(1));

  if(isflat)
    data = zeros([frames length(indexes)]);
  else
    data = zeros([frames ssize(2) length(indexes)]);
  end

  for i=1:frames
    tmpimg = load_data(fname, i);
    for j=1:length(indexes)
      if(isflat)
        data(i, j) = tmpimg(indexes(j));
      else
        data(i,:,j) = tmpimg(indexes(j),:);
        %data(i,:,j) = reshape(tmpimg,[1 ssize(2) 1]);
      end
    end
  end

  %if(length(data)==1)
  %  tmpdata = data(1).cdata;
  %  data = [];
  %  data = tmpdata;
  %end

  return;
end
