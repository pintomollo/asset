function [mymovie] = read_movie(id)
  
  hwait = waitbar(0, 'Extracting the movie','Name','Bio-Formats Library');

  r = loci.formats.ImageReader();
  r.setId(id);

  qt = loci.formats.LegacyQTTools;
  qtjava = qt.canDoQT();
  if (r.getReader().getClass().equals(loci.formats.in.QTReader().getClass()) && qtjava)
    r.getReader().setLegacy(true);
  end

  pixtype = r.getPixelType();
  nbytes = loci.formats.FormatTools.getBytesPerPixel(pixtype);
  imgtype = char(loci.formats.FormatTools.getPixelTypeString(pixtype));

  mymovie = struct('data','');
  %r.setSeries(0);
  w = r.getSizeX();
  h = r.getSizeY();
  numImages = r.getImageCount();
  nRGB = r.getRGBChannelCount();
  
  w = loci.formats.ImageWriter();
  w.setId('test.ome');

  i = 1;
  while i<=numImages
    waitbar(double(i)/(numImages),hwait);
    try
      img = r.openImage(i - 1).getData();
      %pix = r.openBytes(i - 1);
      %size(img)
      %pix = typecast(img,imgtype);
      pix = img.getPixels(0, 0, w, h, []);
      arr = reshape(pix, [w h nRGB]);
      %if nRGB==1
      %  arr = squeeze(arr);
      %end
    catch ME
      disp(ME.message);
      continue;
    end
      
    if i<numImages
      w.saveImage(img.getData(),0,0,0);
    else
      w.saveImage(img.getData(),0,1,1);
    end
    %mymovie.data = store_data(arr, mymovie.data);
    i = i + 1;
  end

  r.close();
  close(hwait);

  return;
end
