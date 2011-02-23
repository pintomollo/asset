function [newfile, policy] = convert_movie(fname, compression, policy)
  %
  % This function convert a movie into a OMETiff movie using the BIO-Formats
  %

  % Import the main classes
  import loci.formats.ImageReader;
  import loci.formats.out.OMETiffWriter;
  import loci.formats.gui.LegacyQTTools;
  import loci.formats.MetadataTools;
  import loci.formats.ChannelMerger;
  import loci.formats.gui.BufferedImageReader;
  import loci.formats.gui.BufferedImageWriter;

  if (nargin < 2)
    compression = '';
    policy = 0;
  elseif (nargin < 3)
    if (ischar(compression))
      policy = 0;
    else
      policy = compression;
      compression = '';
    end
  end

  newfile = '';

  if(isempty(fname) | ~ischar(fname))
    return;
  end

  if (isequal(fname(1),'.'))
    fname = fullfile(pwd, fname(2:end));
  end

  omexmlMeta = MetadataTools.createOMEXMLMetadata();

  % Create the File reader pointing to our file
  r = ImageReader();
  r = ChannelMerger(r);
  r = BufferedImageReader(r);
  r.setMetadataStore(omexmlMeta);

  [tokens,junk]=regexp(fname,'(.+[-_])?([^-_\.]+)(\.[\w\.]+)?','tokens');
  name = tokens{1}{1};
  suffix = tokens{1}{2};
  ext = tokens{1}{3};

  name = [name suffix];

  [slash] = findstr(name,'/');
  if(length(slash)>0)
    slash = slash(end) + 1;
  else
    slash = 1;
  end
  printname = strrep(name(slash:end),'_','\_');

  r.setId(fname);

  % In case it's a QT movie, use the legacy if possible
  qt = LegacyQTTools;
  qtjava = qt.canDoQT();

  if (strncmp(char(r.getFormat()),'QuickTime',9))
      if (qtjava)
        r.getReader().setLegacy(true);
      else
        qt.checkQTLibrary();
      end
  elseif (strncmp(char(r.getFormat()),'OME-TIFF',8))

    newfile = fname;

    try 
      r.close();
    catch ME
      disp(ME.message);
      disp('Ignoring closing Error, Continuing');
    end

    newfile = relativepath(newfile);

    return;
  end
  %r.setId(fname);

  % Create the File writer
  w = OMETiffWriter();
  w = BufferedImageWriter(w);
  w.setMetadataRetrieve(omexmlMeta);

  % Check whether the conversion is feasible
  if(~w.isSupportedType(r.getPixelType()))

    return
  end

  % We create a OMETiff file
  newname = [name '.ome.tif'];

  if(exist(newname,'file'))
    answer = 0;
    if (policy == 0)
      while (answer == 0)
        answer = menu(['The OME-TIFF version of ' printname ' already exists, overwrite it ?'],'Yes', 'Yes to All','No','No to All');
      end
    else
      answer = policy;
    end

    if (mod(answer, 2) == 0)
      policy = answer - 1;
      answer = policy;
    end

    switch answer
      case 1
        delete(newname);
      case 3
        newfile = newname;

        try 
          r.close();
        catch ME
          disp(ME.message);
          disp('Ignoring closing Error, Continuing');
        end

        newfile = relativepath(newfile);

        return;
    end
  end
  w.setId(newname);

  if(isempty(compression))
    % Select the compression
    compress = w.getCompressionTypes();
    [sel,ok] = listdlg('PromptString','Select a Compression',...
                  'SelectionMode','single',...
                  'ListString',char(compress));
    if(ok==0)
      sel = 1;
    end

    compression = compress(sel);
  end
  % Set the compression
  w.setCompression(java.lang.String(compression));

  % A nice status bar
  format = char(r.getReader().getFormat());
  str = [' Converting ' printname ' from [' format '] to [' char(compression) ' OME-TIFF]'];
  hwait = waitbar(0, str,'Name','Bio-Formats Library');
  % Simply copy everything
  numSeries = r.getSeriesCount();
  for s=0:(numSeries-1)
    r.setSeries(s);
    
    numImages = r.getImageCount();

    for i=0:(numImages-1)
      img = r.openImage(i);

      w.saveImage(img, s, (i==(numImages-1)), (i==(numImages-1)) && (s==(numSeries-1)));
      
      % Show our progress
      waitbar(double(s*numSeries + (i+1))/(numSeries*numImages),hwait);
    end
  end

  try 
    r.close();
  catch ME
    disp(ME.message);
    disp('Ignoring Error, Continuing');
  end

  close(hwait);
  
  newfile = newname;
  newfile = relativepath(newfile);

  return;
end
