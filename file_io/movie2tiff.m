function [newfile] = movie2tiff(fname, opts, batch_mode)
% MOVIE2TIFF attempts to convert an input movei recording into a TIFF stack
% using the build-in VideoReader.
%
%   [TIFF] = MOVIE2TIFF(MOV, OPTS) converts MOV into TIFF using VideoReader.
%   If the conversion is not possible, MOVIE2TIFF will display an error
%   message and return MOV instead of TIFF.
%
%   [...] = MOVIE2TIFF(MOV) uses the defaut OPTS structure as returned by
%   get_struct('options');
%
% Wilson lab, University of Otago
% Simon Blanchoud
% 09.06.2015

  % Get the options if need be
  if (nargin < 2)
    opts = get_struct('options');
    batch_mode = false;
  elseif (nargin == 2)
    if (islogical(opts))
      batch_mode = opts;
      opts = get_struct('options');
    else
      batch_mode = false;
    end
  end

  % Split the filename
  [file_path, filename, ext] = fileparts(fname);

  % Open the recording
  try
    video_file = VideoReader(fname);
  catch
    err = lasterror();
    warning(['Error during the video decoding:']);
    print_all(err);

    newfile = fname;
    return;
  end

  % Fancy naming
  printname = strrep(filename,'_','\_');

  % Fancy display
  if (opts.verbosity > 1)
    hwait = waitbar(0, ['Decoding Movie ' printname '...'],'Name','Cell Tracking','Visible','on');
  end

  % Make sure we got the number of frames
  if (~isfield(video_file, 'NumberOfFrames') || ~isfinite(video_file.NumberOfFrames))
    junk = read(video_file, Inf);
  end

  % Create the new file
  newfile = [fname '.tiff'];
  if (exist(newfile, 'file'))

    % We do not accept "empty" answers
    answer = 0;
    if (batch_mode)
      answer = 2;
    end

    while (answer == 0)
      answer = menu(['The TIFF version of ' printname ' already exists, overwrite it ?'],'Yes','No');
    end

    % Act accorindly to the policy
    switch answer

      % Delete the current files (did not dare overwriting it directly)
      case 1
        delete(newfile);

      % Otherwise we can stop here
      case 2
        if (opts.verbosity > 1)
          close(hwait);
        end
        return;
    end
  end

  % Transfer the frames into a TIFF
  prev_img = [];
  for i=1:video_file.NumberOfFrames
    img = read(video_file, i);

    if (isempty(prev_img) || any(img(:) ~= prev_img(:)))
      %try
      %  imwrite(img, newfile, 'TIFF', 'WriteMode', 'append');
      %catch
      %  keyboard
      %end
      newfile = save_stack(newfile, img);
    end

    % Update the loop
    prev_img = img;
    if (opts.verbosity > 1)
      waitbar(i/video_file.NumberOfFrames,hwait);
    end
  end

  if (isstruct(newfile))
    newfile.tiffobj.close();
    newfile = newfile.fname;
  end

  % Close the status bar
  if (opts.verbosity > 1)
    close(hwait);
  end

  return;
end
