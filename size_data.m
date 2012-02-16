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
    if(isstruct(fname) & isfield(fname, 'fname'))
      fname = fname.fname;
    end
    infos = imfinfo(fname);
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

  return;
end
