function [result] = load_data(fname, indexes)
% LOAD_DATA reads image files through imread.
%
%   IMGS = LOAD_DATA(FNAME, INDXS) loads the frames at INDXS from FNAME and returns
%   them as the stack IMGS. IMGS has the same data type as the one used in FNAME.
%   INDXS which are not valid are ignored.
%
%   IMGS = LOAD_DATA(STRUCT, INDXS) utilizes the field 'fname' in STRUCT as FNAME.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 15.05.2014

  % Initialize the output
  result = [];

  if (nargin < 2)
    indexes = [];
  end

  % If this is a structure with proper field, use this file name
  if(isstruct(fname) && isfield(fname, 'fname'))
    fname = fname.fname;
  end

  % Get the stack size
  [nframes, ssize] = size_data(fname);

  % There was something wrong !
  if (nframes < 1)
    return;
  end

  % Remove the invalid indexes
  indexes = indexes(indexes > 0 & indexes <= nframes);

  % Just in case
  if (isempty(indexes))
    return;
  end

  % Read the frames
  result = imread(fname, indexes);
  %% In case we have only one frame to load, we can load it directly
  %if (length(indexes) == 1)
  %  result = imread(fname, indexes);

  %% Otherwise, we need to create the appropriate stack first
  %else
  %  % Load the first frame to know the data type
  %  tmp_img = imread(fname, indexes(1));

  %  % Create the stack
  %  result = zeros([ssize length(indexes)], class(tmp_img));

  %  % Copy the frame
  %  result(:, :, 1) = tmp_img;

  %  % And add the remaining frames
  %  for i = 2:length(indexes)
  %    result(:,:,i) = imread(fname, indexes(i));
  %  end
  %end

  return;
end
