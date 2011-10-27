function [fname] = save_data(fname, images, ntotal)
% STORE_DATA saves into the specified file the provided matrix.
%
%   

  import ome.xml.model.primitives.PositiveInteger;
  import loci.formats.out.OMETiffWriter;
  import loci.formats.MetadataTools;

  if (nargin == 1)
    images = fname;
    fname = [];
    ntotal = size(images, 3);
  elseif (nargin == 2)
    if (isnumeric(fname))
      ntotal = images;
      images = fname;
      fname = [];
    else
      ntotal = size(images, 3);
    end
  end

  if (~isstruct(fname))

    if (isempty(fname))    
      % Get the name of the new file
      fname = absolutepath(get_new_name('tmpmat(\d+)\.ome\.tiff?', 'TmpData'));
    else
      abs_fname = absolutepath(fname);
      mov_fname = absolutepath(fname, 'Movies');
      tmp_fname = absolutepath(fname, 'TmpData');
      if (exist(mov_fname, 'file'))
        fname = mov_fname;
      elseif (exist(tmp_fname, 'file'))
        fname = tmp_fname;
      else  
        fname = abs_fname;
      end
    end

    omexmlMeta = MetadataTools.createOMEXMLMetadata();

    if (~exist(fname))
      [junk, junk, realEndian] = computer;
      bigEndian = ~strncmp(realEndian, 'L', 1);

      omexmlMeta.setImageID('Image:0', 0);
      omexmlMeta.setPixelsID('Pixels:0', 0);
      omexmlMeta.setPixelsSizeX(PositiveInteger(java.lang.Integer(size(images, 2))), 0);
      omexmlMeta.setPixelsSizeY(PositiveInteger(java.lang.Integer(size(images, 1))), 0);
      omexmlMeta.setPixelsSizeZ(PositiveInteger(java.lang.Integer(1)), 0);
      omexmlMeta.setPixelsSizeT(PositiveInteger(java.lang.Integer(ntotal)), 0);
      omexmlMeta.setPixelsSizeC(PositiveInteger(java.lang.Integer(1)), 0);

      omexmlMeta.setPixelsType(ome.xml.model.enums.PixelType.UINT16, 0);
      omexmlMeta.setPixelsBinDataBigEndian(java.lang.Boolean(bigEndian), 0, 0);

      omexmlMeta.setPixelsDimensionOrder(ome.xml.model.enums.DimensionOrder.XYZTC, 0);
      omexmlMeta.setChannelID('Channel:0:0', 0, 0);
      omexmlMeta.setChannelSamplesPerPixel(PositiveInteger(java.lang.Integer(1)), 0, 0);

      indexes = 1:size(images, 3);
      MetadataTools.verifyMinimumPopulated(omexmlMeta, 0)
    else

      error('Cannot add planes at the end of a file (yet?).');

      r = loci.formats.ImageReader();
      r.setMetadataStore(omexmlMeta);
      r.setId(fname);

      nframes = r.getImageCount()
      indexes = nframes + size(images, 3);
      r.close();

      % Update xml ? 
      omexmlMeta.setPixelsSizeT(PositiveInteger(java.lang.Integer(indexes)), 0);

      MetadataTools.verifyMinimumPopulated(omexmlMeta, 0)
    end


    writer = OMETiffWriter();
    writer.setMetadataRetrieve(omexmlMeta);
    writer.setId(fname);
    writer.setCompression(java.lang.String('LZW'))

    %if (~isempty(compress))
    %  writer.setCompression(java.lang.String(compress));
    %end

    offset = 0;
  else
    
    writer = fname.fid;
    offset = fname.offset;
    ntotal = fname.ntotal;
    fname = fname.fname;

      indexes = 1:size(images, 3);
  end

  %keyboard
  [images, img_params] = all2uint16(images);

  for i=1:length(indexes)

    img = images(:,:,i);
    img = img.';
    img = typecast(img(:), 'uint8');

    indx = indexes(i) + offset;

  try
    writer.saveBytes(indx - 1, img);
%  end
catch ME
    print_all(ME);
  end
  end

  if (indx == ntotal)
    writer.close();
    fname = relativepath(fname);
  else
    fname = struct('fname', fname, ...
                   'ntotal', ntotal, ...
                   'fid', writer, ...
                   'offset', indx);
  end

  return;
end
