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

    % Hide the bar in case we are loop several times
    if (opts.verbosity > 1)
      set(hwait, 'Visible','off');
    end

    % Now we extract the corresponding metadata for potential later use
    curdir = pwd;
    cmd_path = which('bfconvert.bat');

    % We need LOCI to do so...
    if (isempty(cmd_path))
      error('ASSET:lociMissing', 'The LOCI command line tools are not present !\nPlease follow the instructions provided by install_cell_tracking');
    end
    [mypath, junk] = fileparts(cmd_path);

    % This can take a while, so inform the user
    hInfo = warndlg('Populating metadata, please wait.', 'Preprocessing movie...');

    % Move to the correct folder
    cd(mypath);

    % And call the LOCI utility to extract the metadata
    if (ispc)
      cmd_name = ['"' fname '"'];
      [res, metadata] = system(['showinf.bat -nopix -nometa -omexml-only ' cmd_name]);

    else
      cmd_name = strrep(fname,' ','\ ');
      [res, metadata] = system(['./showinf -nopix -nometa -omexml-only ' cmd_name]);
    end

    % Delete the information if need be
    if (ishandle(hInfo))
      delete(hInfo);
    end

    % Go back to the original folder
    cd(curdir);

    % Check if an error occured
    if (res ~= 0)
      error('ASSET:Metadata', metadata);
    end

    % Try to identify better metadata
    metadata = find_metadata(fname, metadata);

    % This can take a while, so inform the user
    hInfo = warndlg('Parsing metadata, please wait.', 'Preprocessing movie...');

    % Store the resulting metadata
    [metadata, opts] = parse_metadata(metadata, opts);
    if (length(metadata.channels) == nchannels && nchannels > 1)
      mindx = 0;
      for m = 1:nchannels
        if (~isempty(strfind(fname, metadata.channels{m})))
          mindx = m;
          break
        end
      end

      % Extract the portion corresponding to the current channel
      if (mindx ~= 0)
        metadata.acquisition_time = metadata.acquisition_time(mindx, :);
        metadata.channels = metadata.channels(mindx);
        metadata.exposure_time = metadata.exposure_time(mindx, :);
        metadata.z_position = metadata.z_position(mindx, :);
      end
    end
    myrecording.channels(k).metadata = metadata;

    % Delete the information if need be
    if (ishandle(hInfo))
      delete(hInfo);
    end

    % Store the original file name as we will replace it by the rescaled one
    myrecording.channels(k).file = absolutepath(myrecording.channels(k).fname);

    % Perfom some string formatting for the display
    indx = strfind(myrecording.channels(k).file, filesep);
    if (isempty(indx))
      indx = 1;
    else
      indx = indx(end) + 1;
    end
    if (opts.verbosity > 1)
      waitbar(0, hwait, ['Preprocessing Movie ' strrep(myrecording.channels(k).file(indx:end),'_','\_')]);
      set(hwait, 'Visible', 'on');
    end

    % Get the absolute file name
    fname = myrecording.channels(k).file;

    % Get the name of the new file
    tmp_fname = absolutepath(get_new_name('tmpmat(\d+)\.ome\.tiff?', 'TmpData'));

    % Get the number of frames
    nframes = size_data(fname);

    % Temporary parameters about the type of data contained in the reader
    img_params = [];

    % Loop over the frames
    for i=1:nframes
      % Convert the image into UINT16
      [img, img_params] = all2uint16(load_data(fname, i), img_params);

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
      if(minimg < myrecording.channels(k).min)
        myrecording.channels(k).min = minimg;
      end
      if(maximg > myrecording.channels(k).max)
        myrecording.channels(k).max = maximg;
      end

      % Save the image in the temporary file
      tmp_fname = save_stack(tmp_fname, img);

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
        tmp_fname = save_stack(tmp_fname, img);

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

function metadata = find_metadata(filename, metadata)
% This function tries to identify more suitable metadata. For now
% on, the following metadata are supported:
%   - Leica Application Suite ".las"
%   - Leica Application Suite ".cal.xml"
%   - uManager "metadata.txt"
%   - files manually placed in the "Metadata" folder and named
%     as the recording

  % Get the folder in which the file is contained
  [file_path, file_name, file_ext] = fileparts(filename);

  % Check if the file exists
  if (exist(fullfile(file_path, '.las'), 'file'))

    % Load it !
    metadata = fileread(fullfile(file_path, '.las'));

  % The other LAS
  elseif (exist(fullfile(file_path, '.Metadata'), 'dir'))
    if (exist(fullfile(file_path, '.Metadata', [file_name file_ext '.cal.xml']), 'file'))
      metadata = fileread(fullfile(file_path, '.Metadata', [file_name file_ext '.cal.xml']));
    else
      indx = strfind(file_name, '_');
      if (~isempty(indx))
        file = regexpdir(fullfile(file_path, '.Metadata'), [file_name(1:indx(end)-1) '(\..*)?\.cal\.xml' ]);
        if (length(file)==1)
          metadata = fileread(file{1});
        end
      end
    end

  % For uManager
  elseif (exist(fullfile(file_path, 'metadata.txt'), 'file'))
    metadata = fileread(fullfile(file_path, 'metadata.txt'));

  % For manually edited files
  elseif (exist(fullfile(pwd, 'Metadata', [filename '.txt']), 'file'))
    metadata = fileread(fullfile(pwd, 'Metadata', [filename '.txt']));
  end

  return;
end
