function tmp_fname = tmp2tiff(fname, template)

  [imgsize, nframes] = size_tmp(fname);
   
  %omexmlMeta = loci.formats.MetadataTools.createOMEXMLMetadata();
  %reader = loci.formats.ImageReader();
  %reader = loci.formats.ChannelMerger(reader);
  %reader.setMetadataStore(omexmlMeta);
  %reader.setId(template);
  %reader.close();

  %omexmlMeta.setPixelsSizeT(ome.xml.model.primitives.PositiveInteger(java.lang.Integer(nframes)), 0);

  tmp_fname = get_temp_name();

  %writer = loci.formats.out.OMETiffWriter();
  %writer.setMetadataRetrieve(omexmlMeta);
  %writer.setId(tmp_fname);
  %writer.setCompression(java.lang.String('LZW'));
  %writer.setWriteSequentially(true);

  for i=1:nframes
    img = load_tmp(fname, i);
    %store_data(writer, img, i);
    imwrite(img, tmp_fname, 'TIFF', 'WriteMode', 'append');
  end

  tmp_fname = relativepath(tmp_fname);
  %writer.close();

  return;
end
