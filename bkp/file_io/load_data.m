function [result] = load_data(fid, indexes)
  % LOAD_DATA reads the requested frames from a file using
  %   imread.

  if(isstruct(fid) & isfield(fid, 'fname'))
    fid = fid.fname;
  end

  [nframes, ssize] = size_data(fid);

  indexes = indexes(indexes > 0 & indexes <= nframes);

  if (length(indexes) == 1)
    result = imread(fid, indexes);
  else
    % If RGB we can have some problems 
    %result = imread(fid, indexes);

    tmp_img = imread(fid, indexes(1));
    nchannels = size(tmp_img, 3);
    result = zeros([ssize nchannels length(indexes)], class(tmp_img));
    result(:, :, :, 1) = tmp_img;

    for i = 2:length(indexes)
      result(:,:,:,i) = imread(fid, indexes(i));
    end
    result = squeeze(result);
  end

  return;
end
