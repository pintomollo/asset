function [nframes, ssize, pixelType, r] = size_data(fname)
  
  if (~isjava(fname))
    if(isstruct(fname) & isfield(fname, 'fname'))
      fname = fname.fname;
    end

    %omexmlMeta = loci.formats.MetadataTools.createOMEXMLMetadata();

    r = loci.formats.ImageReader();
    %r = loci.formats.ChannelMerger(r);
    %r = loci.formats.ChannelSeparator(r);
    %r.setMetadataStore(omexmlMeta);
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

  return;

%  h = fopen(fname,'r+');
%
%  frewind(h);
%  ssize = fread(h,3,'double')';
%  datasize = prod(ssize);
%  doublesize = ftell(h) / 3;
%
%  fseek(h,0,'eof');
%  frames = ftell(h);
%  frames = ((frames / doublesize) - 3) / datasize;
%
%  frame_size = datasize*doublesize;
%
%  if(nargout>3)
%    fseek(h,doublesize*3,'bof');
%  else
%    fclose(h);
%  end
%
%  return;
end
