function [mymovie] = rescale_movie(mymovie, overwrite)
  
  fields = fieldnames(mymovie);
  [junk, junk, realEndian] = computer;
  bigEndian = ~strncmp(realEndian, 'L', 1);

  hwait = waitbar(0,'','Name','CellCoord Info');
  for f = 1:length(fields)
    field = fields{f};

    for k=1:length(mymovie.(field))
      mymovie.(field)(k).file = mymovie.(field)(k).fname;
      indx = findstr(mymovie.(field)(k).file, '/');
      if (isempty(indx))
        indx = 0;
      else
        indx = indx(end) + 1;
      end
      waitbar(0, hwait, ['Preprocessing Movie ' strrep(mymovie.(field)(k).file(indx:end),'_','\_')]);

      tmp_fname = get_new_name('tmpmat(\d+)\.ome\.tiff?', 'TmpData');
      %tmp_fname = get_temp_name();
      %mymovie.(field)(k).fname

      %keyboard

      omexmlMeta = loci.formats.MetadataTools.createOMEXMLMetadata();
      reader = loci.formats.ImageReader();
      reader = loci.formats.ChannelMerger(reader);
      reader.setMetadataStore(omexmlMeta);
      reader.setId(mymovie.(field)(k).file);
      %reader.setSeries(0);
      
      omexmlMeta.setPixelsType(ome.xml.model.enums.PixelType.DOUBLE, 0);
      omexmlMeta.setPixelsBinDataBigEndian(java.lang.Boolean(bigEndian), 0, 0);

      writer = loci.formats.out.OMETiffWriter();
      writer.setMetadataRetrieve(omexmlMeta);
      writer.setId(tmp_fname);
      writer.setCompression(java.lang.String(mymovie.(field)(k).compression));
      writer.setWriteSequentially(true);

      nframes = reader.getImageCount();
      %new_field = mymovie.(field)(k);
      %new_field.fname = '';

      for i=1:nframes
        img = double(load_data(reader, i));


        if (i==1)
          orig = img;
        end

        if (mymovie.(field)(k).hot_pixels)
          img = imhotpixels(img);
        end
        if (mymovie.(field)(k).detrend)
          img = imdetrend(img);
        end


    %img = imnorm(img);
    %img = img.';
    %img = typecast(img(:), 'uint8');
    %writer.saveBytes(i - 1, img);

        minimg = min(img(:));
        maximg = max(img(:));

        if(minimg < mymovie.(field)(k).min)
          mymovie.(field)(k).min = minimg;
        end
        if(maximg > mymovie.(field)(k).max)
          mymovie.(field)(k).max = maximg;
        end
        
        writer = store_data(writer, img, i);

        if (i==1)
          worked = img;
        end

        waitbar(i/(2*nframes),hwait);
      end

      reader.close();
      writer.close();

      mymovie.(field)(k).fname = get_new_name('tmpmat(\d+)\.ome\.tiff?', 'TmpData');
      %mymovie.(field)(k).fname = get_temp_name();

      %clear reader writer;

      %omexmlMeta = loci.formats.MetadataTools.createOMEXMLMetadata();
      %reader = loci.formats.ImageReader();
      %reader = loci.formats.ChannelMerger(reader);
      %reader.setMetadataStore(omexmlMeta);
      reader.setId(tmp_fname);
      
      %omexmlMeta.setPixelsType(ome.xml.model.enums.PixelType.DOUBLE, 0);

      %writer = loci.formats.out.OMETiffWriter();
      %writer.setMetadataRetrieve(omexmlMeta);
      writer.setId(mymovie.(field)(k).fname);
      writer.setCompression(java.lang.String(mymovie.(field)(k).compression));
      writer.setWriteSequentially(true);



      %reader.setMetadataStore(omexmlMeta);
      %reader.setId(tmp_fname);

      %writer.setMetadataRetrieve(omexmlMeta);
      %writer.setId(mymovie.(field)(k).fname);
      %writer.setCompression(java.lang.String(mymovie.(field)(k).compression));
      %writer.setWriteSequentially(true);

%      keyboard
      %try
      for i=1:nframes

        img = load_data(reader, i);
        %reader.close();

        img = imnorm(img, mymovie.(field)(k).min, mymovie.(field)(k).max);

        %writer.setId(mymovie.(field)(k).fname);
        %writer.setCompression(java.lang.String(mymovie.(field)(k).compression));
        store_data(writer, img, i);
        %writer.close();

        waitbar(0.5 + i/(2*nframes),hwait);
      end
      %catch ME
      %  beep;keyboard
      %end

      reader.close();
      writer.close();

%      keyboard

      delete(tmp_fname);
    end
  end

  close(hwait);

  clear img;

  return;
end


