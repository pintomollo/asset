function [myrecording, opts] = preprocess_movie(myrecording, opts)
% PREPROCESS_MOVIE converts the OME-TIFF recordings contained in a tracking structure
% into properly filtered (as defined by the structure, see inspect_channels.m) UINT16
% files.
%
%   [MYRECORDING] = PREPROCESS_MOVIE(MYRECORDING, OPTS) rescales all the recordings used
%   in the tracking experiment using OPTS.
%
%   [MYRECORDING, OPTS] = PREPROCESS_MOVIE(...) returns in addition OPTS.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 18.05.2014

  % Initialize some computer-specific variables required for the conversion
  maxuint = intmax('uint16');
  minuint = intmin('uint16');

  % A nice status-bar if possible
  if (opts.verbosity > 1)
    hwait = waitbar(0,'','Name','Cell Tracking','Visible','off');
  end

  % Get the number of channels to parse
  nchannels = length(myrecording.channels);

  % Loop over all of them
  for k = 1:nchannels

    % We need the absolute path for Java to work properly
    fname = absolutepath(myrecording.channels(k).fname);

    % Hide the bar in case we are looping several times
    if (opts.verbosity > 1)
      set(hwait, 'Visible','off');
    end

    % Split the metadata between the respective channels
    metadata = myrecording.channels(k).metadata;
    if (length(metadata.channels) == nchannels && nchannels > 1)
      mindx = 0;
      for m = 1:nchannels
        if (~isempty(strfind(fname, metadata.channels{m})) || ~isempty(strfind(fname, strrep(metadata.channels{m}, ' ', '-'))))
          mindx = m;
          break
        end
      end

      % Extract the portion corresponding to the current channel
      if (mindx ~= 0)
        [c,p,f] = size(metadata.acquisition_time);
        findx = repmat([1:f], p, 1);
        pindx = repmat([1:p].', 1, f);

        metadata.acquisition_time = metadata.acquisition_time(mindx, :, :);
        metadata.channels = metadata.channels(mindx);
        metadata.exposure_time = metadata.exposure_time(mindx, :, :);
        metadata.z_position = metadata.z_position(mindx, :, :);
        metadata.files = metadata.files(mindx, :, :);

        metadata.acquisition_time = metadata.acquisition_time(:).';
        metadata.exposure_time = metadata.exposure_time(:).';
        metadata.z_position = metadata.z_position(:).';
        metadata.files = metadata.files(:).';

        metadata.frame_index = findx(:).';
        metadata.plane_index = pindx(:).';
      end
    end
    myrecording.channels(k).metadata = metadata;

    % Store the original file name as we will replace it by the rescaled one
    myrecording.channels(k).file = fname;

    % Perfom some string formatting for the display
    [junk, tmp_name, junk] = fileparts(fname);
    [junk, tmp_name, junk] = fileparts(tmp_name);
    if (opts.verbosity > 1)
      waitbar(0, hwait, ['Preprocessing Movie ' strrep(tmp_name,'_','\_')]);
      set(hwait, 'Visible', 'on');
    end

    % Get the name of the new file
    tmp_fname = absolutepath(get_new_name('tmpmat(\d+)\.ome\.tiff?', 'TmpData'));

    % Get the number of frames
    nframes = size_data(fname);

    % Temporary parameters about the type of data contained in the reader
    img_params = [];

    % Loop over the frames
    for i=1:nframes
      % Convert the image into UINT16
      [img, img_params] = scaled_cast(load_data(fname, i), img_params);

      % Perform the required filtering
      if (myrecording.channels(k).detrend)
        img = imdetrend(img, opts.filtering.detrend_meshpoints);
      end
      if (myrecording.channels(k).cosmics)
        img = imcosmics(img, opts.filtering.cosmic_rays_window_size, opts.filtering.cosmic_rays_threshold);
      end
      if (myrecording.channels(k).hot_pixels)
        img = imhotpixels(img, opts.filtering.hot_pixels_threshold);
      end

      % Get the current range of values
      minimg = min(img(:));
      maximg = max(img(:));

      % We'll store the biggest range, to rescale it afterwards
      if (minimg < myrecording.channels(k).min)
        myrecording.channels(k).min = minimg;
      end
      if (maximg > myrecording.channels(k).max)
        myrecording.channels(k).max = maximg;
      end

      % Save the image in the temporary file
      tmp_fname = save_data(tmp_fname, img);

      % Update the progress bar if needed
      if (opts.verbosity > 1)
        if (myrecording.channels(k).normalize)
          waitbar(i/(2*nframes),hwait);
        else
          waitbar(i/nframes,hwait);
        end
      end
    end

    if (isstruct(tmp_fname))
      tmp_fname.tiffobj.close();
      tmp_fname = tmp_fname.fname;
    end

    % Rescale if required by the user
    if (myrecording.channels(k).normalize)
      % Get a third file to write into
      fname = tmp_fname;
      tmp_fname = absolutepath(get_new_name('tmpmat(\d+)\.ome\.tiff?', 'TmpData'));
      myrecording.channels(k).fname = tmp_fname;

      % Loop again over the frames
      for i=1:nframes

        % Load and rescale using the previously measured range
        img = load_data(fname, i);
        img = imnorm(img, myrecording.channels(k).min, myrecording.channels(k).max, '', 0, maxuint);

        % Save the file
        tmp_fname = save_data(tmp_fname, img);

        % Update the progress bar
        if (opts.verbosity > 1)
          waitbar(0.5 + i/(2*nframes),hwait);
        end
      end

      if (isstruct(tmp_fname))
        tmp_fname.tiffobj.close();
        tmp_fname = tmp_fname.fname;
      end

      % Delete the intermidary file (i.e. the filtered one)
      delete(fname);
    else
      myrecording.channels(k).fname = tmp_fname;
    end

    % Get everything in relative paths
    myrecording.channels(k).file = relativepath(myrecording.channels(k).file);
    myrecording.channels(k).fname = relativepath(myrecording.channels(k).fname);
  end

  % Close the status bar
  if (opts.verbosity > 1)
    close(hwait);
  end

  return;
end
