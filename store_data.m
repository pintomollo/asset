function [fname] = store_data(fname, images, indexes)
% STORE_DATA saves into the specified file the provided matrix.
%
%   

  if (~isjava(fname))

    compress = '';
    if (isstruct(fname) & isfield(fname, 'fname'))

      if (isfield(fname, 'compression'))
        compress = fname.compression;
      end

      fname = fname.fname;
    end
    omexmlMeta = loci.formats.MetadataTools.createOMEXMLMetadata();

    fname = absolutepath(fname);

    r = loci.formats.ImageReader();
    r.setMetadataStore(omexmlMeta);
    r.setId(fname);

    if (nargin < 3)
      nframes = r.getImageCount();
      indexes = nframes + 1;

      % Update xml ? 
      omexmlMeta.setPixelsSizeT(ome.xml.model.primitives.PositiveInteger(java.lang.Integer(indexes)), 0);
    end

    r.close();

    writer = loci.formats.ImageWriter();
    writer.setMetadataRetrieve(omexmlMeta);
    writer.setId(fname);

    if (~isempty(compress))
      writer.setCompression(java.lang.String(compress));
    end
  else
    writer = fname;
  end

  %keyboard

  for i=1:length(indexes)

    img = images(:,:,i);
    img = img.';
    img = typecast(img(:), 'uint8');

    writer.saveBytes(indexes(i) - 1, img);
  end

  if (~isjava(fname))
    writer.close();
  end

  return;
end
