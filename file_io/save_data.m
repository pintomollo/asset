function done = save_data(fname, imgs)
% SAVE_DATA stores images into the provided filename as stack TIFF files using
% imwrite.
%
%   DONE = SAVE_DATA(FNAME, IMG) stores IMG in FNAME as a stack TIFF file. If FNAME
%   does not exist, it creates it. If it does, IMG is appended at the end of it.
%   DONE is true if the saving worked properly.
%
%   DONE = SAVE_DATA(FNAME, STACK) stores the whole STACK in FNAME.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 15.05.2014

  done = false;

  % If no filename is provided, we have a problem
  if (isempty(fname))
    warning('CAST:save_data', 'No file name was provided.');
    return;

  % Or, potentially is it's a structure
  elseif(isstruct(fname))

    % We're looking only for one field with the name 'fname', everyhting else if wrong
    if (~isfield(fname, 'fname') || isempty(fname.fname))
      warning('CAST:save_data', 'No file name was found in the provided structure.');
      return;

    % We'll try this one then
    else
      fname = fname.fname;
    end

  % If we have a string, try to open the file
  elseif (~ischar(fname))
    warning('CAST:save_data', 'Unable to extract a filename from an "%s" object.', class(fname));
    return
  end

  % Check that we can actually save what we got
  if (isempty(imgs) || ~isnumeric(imgs))
    warning('CAST:save_data', 'No valid image was provided.');
    return
  end

  imwrite(imgs, fname, 'TIFF', 'WriteMode', 'append');
  done = true;
  %% Loop over the planes to save all of them
  %nplanes = size(imgs, 3);
  %for i=1:nplanes

  %  % The current image
  %  img = imgs(:,:,i);

  %  % Now here is something incredible, this short "waiting" hack was necessary as,
  %  % on a high-performance PC, Matlab was try to save too fast in the file and
  %  % got messed-up with the handlers. Not sure it is still necessary though.
  %  waiting_time = 0;
  %  while (~done)

  %    % Try to save the plane in the file
  %    try
  %      imwrite(img, fname, 'TIFF', 'WriteMode', 'append');

  %      % Done then !
  %      done = true;

  %    % Here we have a problem, then we try again in a few time, might be better !
  %    catch ME

  %      % In case we have not waited so much, we try again
  %      if (waiting_time < 20)
  %        nsecs = rand(1);
  %        waiting_time = waiting_time + nsecs;
  %        pause(nsecs);

  %      % Otherwise, given up
  %      else
  %        rethrow(ME)
  %      end
  %    end
  %  end
  %end

  return
end
