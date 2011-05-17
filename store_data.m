function [fid] = store_data(fid, images, indexes)

  if (~isjava(fid))

    compress = 'LZW';
    if (isstruct(fid) & isfield(fid, 'fname'))

      if (isfield(fid, 'compression'))
        compress = fid.compression;
      end

      fid = fid.fname;
    end
    omexmlMeta = loci.formats.MetadataTools.createOMEXMLMetadata();

    fid = absolutepath(fid);

    r = loci.formats.ImageReader();
    r = loci.formats.ChannelMerger(r);
    r.setMetadataStore(omexmlMeta);
    r.setId(fid);

    if (nargin < 3)
      nframes = r.getImageCount();
      indexes = nframes + 1;

      % Update xml ? 
      omexmlMeta.setPixelsSizeT(ome.xml.model.primitives.PositiveInteger(java.lang.Integer(indexes)), 0);
    end

    r.close();

    writer = loci.formats.ImageWriter();
    writer.setMetadataRetrieve(omexmlMeta);
    writer.setId(fid);

    writer.setCompression(java.lang.String(compress));
  else
    writer = fid;
  end

  %keyboard

  for i=1:length(indexes)

    img = images(:,:,i);
    img = img.';
    img = typecast(img(:), 'uint8');

    writer.saveBytes(indexes(i) - 1, img);
  end

  if (~isjava(fid))
    writer.close();
  end

  return;
end
