function [nframes, ssize, pixelType, r] = size_data(fname)
  
  if (nargout == 4 | isjava(fname))
    if (~isjava(fname))
      if(isstruct(fname) & isfield(fname, 'fname'))
        fname = fname.fname;
      end

      fname = absolutepath(fname);

      r = loci.formats.ImageReader();
      r.setId(fname);
    else
      r = fname;
    end

    height = r.getSizeX();
    width = r.getSizeY();
    ssize = [width height];
    nframes = r.getImageCount();

    if (nargout > 2)
      pixelType = r.getPixelType();
    end
    if (nargout < 4 && ~isjava(fname))
      r.close();
    end
  else
    %Substitute imfinfos to handle natively missing files
    if(isstruct(fname) & isfield(fname, 'fname'))
      fid = fopen(fname.fname, 'r');
    else
      %infos = imfinfo(fname);
      fid = fopen(fname, 'r');
    end

    if (fid == -1)
      if (isstruct(fname))
        if (isfield(fname, 'eggshell'))
          nframes = length(fname.eggshell);
        elseif (isfield(fname, 'cortex'))
          nframes = length(fname.cortex);
        elseif (isfield(fname, 'centers'))
          nframes = size(fname.centers, 2);
        else
          error('MATLAB:imfinfo:fileOpen', ...
              'Unable to find suitable fields to get the size of the recording.');
        end
      else
        error('MATLAB:imfinfo:fileOpen', ...
            'Unable to open file "%s" for reading.', fname);
      end

      ssize = NaN(1,2);
      pixelType = 'unknown';
    else
      filename = fopen(fid);  % Get the full pathname if not in pwd.
      fclose(fid);

      % Determine filetype from file.
      [format, fmt_s] = imftype(filename);
      infos = feval(fmt_s.info, filename);

      nframes = length(infos);
      ssize = [infos(1).Height infos(1).Width];
      if (nargout > 2)
        switch lower(infos(1).SampleFormat(1))
          case 'u'
            pixelType = ['uint' num2str(infos(1).BitDepth)];
          case 'i'
            pixelType = ['int' num2str(infos(1).BitDepth)];
          otherwise
            pixelType = 'double';
        end
      end
    end
  end

  return;
end
