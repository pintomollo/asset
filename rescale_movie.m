function [mymovie] = rescale_movie(mymovie, opts)
% RESCALE_MOVIE converts an OME-TIFF recording into a UINT16 one, rescaling its values
% to use the available range of values. In addition it performs the required filtering
% (i.e. hot-pixels removal and detrending).
%
%   [MYMOVIE] = RESCALE_MOVIE(MYMOVIE, OPTS) rescales all the recordings used in the
%   experiement (i.e. each channel of the recording) using the options OPTS.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 20.06.2011

  % Initialize some computer-specific variables required for the conversion
  %[junk, junk, realEndian] = computer;
  %bigEndian = ~strncmp(realEndian, 'L', 1);
  maxuint = intmax('uint16');

  % Import the Java classes
  %import loci.formats.ImageReader;
  %import loci.formats.out.OMETiffWriter;
  %import loci.formats.MetadataTools;
  %import loci.formats.ChannelMerger;

  % A nice status-bar if possible
  if (opts.verbosity > 1)
    hwait = waitbar(0,'','Name','CellCoord Info');
  end

  % Get all the fields of the experiement as this might change based on the data,
  % and loop over them
  fields = fieldnames(mymovie);
  for f = 1:length(fields)
    field = fields{f};

    % If the current field does not contain a file, skip it
    if (~isfield(mymovie.(field), 'fname'))
      continue;
    end
  
    % Fields can be arrays of structures, so loop over them
    for k = 1:length(mymovie.(field))
      
      % Stire the original file name as we will replace it by the rescaled one
      mymovie.(field)(k).file = mymovie.(field)(k).fname;

      % Perfom some string formatting for the display
      indx = strfind(mymovie.(field)(k).file, filesep);
      if (isempty(indx))
        indx = 1;
      else
        indx = indx(end) + 1;
      end
      if (opts.verbosity > 1)
        waitbar(0, hwait, ['Preprocessing Movie ' strrep(mymovie.(field)(k).file(indx:end),'_','\_')]);
      end

      fname = absolutepath(mymovie.(field)(k).file);

      % Get the name of the new file
      tmp_fname = absolutepath(get_new_name('tmpmat(\d+)\.ome\.tiff?', 'TmpData'));
  
      % Create the metadata container
      %omexmlMeta = MetadataTools.createOMEXMLMetadata();

      % Create the reader
      %reader = ImageReader();
      % If required, try to merge separated files
      %if (opts.merge_input_files)
      %  r = ChannelMerger(reader);
      %end
      % Initialize it
      %reader.setMetadataStore(omexmlMeta);
      %reader.setId(absolutepath(mymovie.(field)(k).file));
      
      % Set the writer rescaled-specific parameters
      %omexmlMeta.setPixelsType(ome.xml.model.enums.PixelType.UINT16, 0);
      %omexmlMeta.setPixelsBinDataBigEndian(java.lang.Boolean(bigEndian), 0, 0);

      % Create the writer
      %writer = OMETiffWriter();
      %writer.setMetadataRetrieve(omexmlMeta);
      %writer.setId(tmp_fname);
      %writer.setCompression(java.lang.String(mymovie.(field)(k).compression));
      %writer.setWriteSequentially(true);

      % Get the number of frames
      %nframes = reader.getImageCount();
      nframes = size_data(fname);

      % Temporary parameters about the type of data contained in the reader
      img_params = [];

      % Loop over the frames
      for i=1:nframes
        % Convert the image into UINT16
        [img, img_params] = all2uint16(load_data(fname, i), img_params);

        % Perform the required filtering
        if (mymovie.(field)(k).hot_pixels)
          img = imhotpixels(img);
        end
        if (mymovie.(field)(k).detrend)
          img = imdetrend(img);
        end

        % Get the current range of values
        minimg = min(img(:));
        maximg = max(img(:));

        % We'll store the biggest range, to rescale it afterwards
        if(minimg < mymovie.(field)(k).min)
          mymovie.(field)(k).min = minimg;
        end
        if(maximg > mymovie.(field)(k).max)
          mymovie.(field)(k).max = maximg;
        end
        
        % Write the filtered data, we'll rescale them in a second pass as we do not
        % know the range yet
        %writer = store_data(writer, img, i);
        done = false;
        waiting_time = 0;
        while (~done)
          try
            imwrite(img, tmp_fname, 'TIFF', 'WriteMode', 'append');
            done = true;
          catch ME
            if (waiting_time < 20)
              nsecs = rand(1);
              waiting_time = waiting_time + nsecs;
              pause(nsecs);
            else
              rethrow(ME)
            end
          end
        end

        % Update the progress bar if needed
        if (opts.verbosity > 1)
          waitbar(i/(2*nframes),hwait);
        end
      end

      % Close both handlers
      %reader.close();
      %writer.close();

      % Get a third file to write into
      fname = tmp_fname;
      tmp_fname = absolutepath(get_new_name('tmpmat(\d+)\.ome\.tiff?', 'TmpData'));
      mymovie.(field)(k).fname = tmp_fname;

      % Reade from the new filtered file
      %reader.setId(tmp_fname);

      % Create a new writer
      %writer = OMETiffWriter();
      %writer.setMetadataRetrieve(omexmlMeta);
      %writer.setId(absolutepath(mymovie.(field)(k).fname));
      %writer.setCompression(java.lang.String(mymovie.(field)(k).compression));
      %writer.setWriteSequentially(true);

      % Loop again over the frames
      for i=1:nframes

        % Load and rescale using the previously measured range
        try
        img = load_data(fname, i);
        catch
          beep;keyboard
        end
        img = imnorm(img, mymovie.(field)(k).min, mymovie.(field)(k).max, '', 0, maxuint);

        % And save the final image
        done = false;
        waiting_time = 0;
        while (~done)
          try
            imwrite(img, tmp_fname, 'TIFF', 'WriteMode', 'append');
            done = true;
          catch ME
            if (waiting_time < 20)
              nsecs = rand(1);
              waiting_time = waiting_time + nsecs;
              pause(nsecs);
            else
              rethrow(ME)
            end
          end
        end
        %store_data(writer, img, i);

        % Update the progress bar
        if (opts.verbosity > 1)
          waitbar(0.5 + i/(2*nframes),hwait);
        end
      end

      % Close both handlers
      %reader.close();
      %writer.close();

      % Delete the intermidary file (i.e. the filtered one)
      delete(fname);
    end
  end

  % Close the status bar
  if (opts.verbosity > 1)
    close(hwait);
  end

  % Clean the memory
  clear img;

  return;
end
