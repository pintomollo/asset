function export_movie(mymovie, opts)

  if (nargin < 2)
    opts = struct('compression', '', ...
                  'max_export', 1, ...
                  'crop_export', false, ...
                  'parse_export', 'normal');

  end

  import loci.common.DataTools;
  import loci.formats.gui.AWTImageTools;
  import loci.formats.gui.BufferedImageWriter;
  import loci.formats.out.OMETiffWriter;
  import loci.formats.MetadataTools;

  colorModel = AWTImageTools.makeColorModel(1, java.awt.image.DataBuffer.TYPE_BYTE);
  policy = 0;
  do_init = true;

  fields = fieldnames(mymovie);
  for f = 1:length(fields)
    field = fields{f};

    for k = 1:length(mymovie.(field))

      if (~isfield(mymovie.(field)(k), 'file'))
        continue;
      else
        crop_export = opts.crop_export;

        if (crop_export)
          if (isfield(mymovie.(field)(k), 'centers'))
            lookup_field = field;
            lookup_indx = k;
          elseif (isfield(mymovie, 'markers') && isfield(mymovie.markers, 'centers'))
            lookup_field = 'markers';
            lookup_indx = 1;
          else
            crop_export = false;
            warning(['No segmentation data allows to crop the file ' mymovie.(field)(k).file]);
          end
        end


        [tokens,junk]=regexp(mymovie.(field)(k).file,'(.+[-_])?([^-_\.]+)(\.[\w\.]+)?','tokens');
        name = tokens{1}{1};
        suffix = tokens{1}{2};
        ext = tokens{1}{3};

        if (~strcmp(ext,'.ome.tif'))
          ext = '.ome.tif';
        end

        if (crop_export)
          diff_name = 'CROP';
        else
          diff_name = 'RECOS';
        end

        if (~isempty(name))
          new_file = [name diff_name '-' suffix ext];
        else
          new_file = [suffix '-' diff_name ext];
        end

        [slash] = findstr(new_file,'/');
        if(length(slash)>0)
          slash = slash(end) + 1;
        else
          slash = 1;
        end
        printname = strrep(new_file(slash:end),'_','\_');

        if(exist(new_file,'file'))
          answer = 0;
          if (policy == 0)
            while (answer == 0)
              answer = menu([printname ' already exists, overwrite it ?'],'Yes','Yes to All','No','No to All');
            end
            drawnow;
          else
            answer = policy;
          end

          if (mod(answer, 2) == 0)
            policy = answer - 1;
            answer = policy;
          end

          switch answer
            case 1
              delete(new_file);
            case 3
              try 
                r.close();
              catch ME
                disp(ME.message);
                disp('Ignoring closing Error, Continuing');
              end

              continue;
          end
        end
      end

      if (crop_export)

        if (length(opts.crop_size) == 1)
          crop_size = round(2 * opts.crop_size * max(mymovie.(lookup_field)(lookup_indx).axes_length([2 1], :), [], 2));
        else
          crop_size = opts.crop_size;
          if (any(crop_size(:) < max(mymovie.(lookup_field)(lookup_indx).axes_length([2 1], :), [], 2)))
            warning('Crop size is smaller than some estimated embryos''s size')
          end
        end
      end
      %import loci.formats.gui.ImageViewer;

      if (do_init)
        [imgsize nframes] = size_data(mymovie.(field)(k));

        switch (opts.parse_export)
          case 'normal'
            frames = 1:nframes; 
          case 'random'
            frames = randperm(nframes);      
        end

        if (opts.max_export > 1)
          max_frames = round(opts.max_export);
        else
          max_frames = round(nframes*opts.max_export);
        end

        if max_frames > nframes 
          max_frames = nframes;
        elseif max_frames < 1
          max_frames = 1;
        end

        frames = frames(1:max_frames) - 1;
        do_init = ~do_init;
      end

      % We convert an image to get the corresponding values for the metadata
      if (isfield(mymovie.(field)(k), 'metadata'))
        omexmlMetaOrig = MetadataTools.createOMEXMLMetadata(mymovie.(field)(k).metadata);

        new_xml = char(mymovie.(field)(k).metadata);
        [starts, ends] = regexp(new_xml,'<Plane.*?>.*?</\Plane>');
        new_xml = new_xml([1:ends(max_frames) (ends(end)+1):length(new_xml)]);

        omexmlMeta = MetadataTools.createOMEXMLMetadata(java.lang.String(new_xml));
        omexmlMeta.setPixelsPixelType(java.lang.String('uint8'), 0, 0);

        has_metadata = true;
      else
        % could extend to export whole movie to a single file
        channels = 1;

        omexmlMeta = MetadataTools.createOMEXMLMetadata();
        omexmlMeta.createRoot();

        omexmlMeta.setPixelsSizeX(java.lang.Integer(imgsize(2)), 0, 0);
        omexmlMeta.setPixelsSizeY(java.lang.Integer(imgsize(1)), 0, 0);
        omexmlMeta.setPixelsSizeZ(java.lang.Integer(1), 0, 0);
        omexmlMeta.setPixelsSizeC(java.lang.Integer(channels), 0, 0);
        omexmlMeta.setPixelsPixelType(java.lang.String('uint8'), 0, 0);
        omexmlMeta.setPixelsBigEndian(java.lang.Boolean(false), 0, 0);
        omexmlMeta.setPixelsDimensionOrder(java.lang.String('XYCZT'), 0, 0);
        omexmlMeta.setLogicalChannelSamplesPerPixel(java.lang.Integer(channels), 0, 0);

        has_metadata = false;
      end

      if (crop_export)
        omexmlMeta.setPixelsSizeX(java.lang.Integer(crop_size(2)), 0, 0);
        omexmlMeta.setPixelsSizeY(java.lang.Integer(crop_size(1)), 0, 0);
      end

      omexmlMeta.setPixelsSizeT(java.lang.Integer(max_frames), 0, 0);
      for i = 1:max_frames
        nframe = frames(i);
        if (has_metadata)
          omexmlMeta.setPlaneTimingDeltaT(omexmlMetaOrig.getPlaneTimingDeltaT(0, 0, nframe), 0, 0, i - 1);
        else
          omexmlMeta.setPlaneTimingDeltaT(java.lang.Double(nframes * 10), 0, 0, i - 1);
        end
      end

      if (isequal(new_file(1),'.'))
        new_file = fullfile(pwd, new_file(2:end));
      end

      w = OMETiffWriter();
      w = BufferedImageWriter(w);
      w.setMetadataRetrieve(omexmlMeta);
      w.setId(new_file);
      w.setInterleaved(false);

      if(isempty(opts.compression))
        % Select the compression
        compression = w.getCompressionTypes();
        [sel,ok] = listdlg('PromptString','Select a Compression',...
                      'SelectionMode','single',...
                      'ListString',char(compression));
        if(ok==0)
          sel = 1;
        end

        opts.compression = compression(sel);
      end
      % Set the compression
      w.setCompression(java.lang.String(opts.compression));

      %viewer = ImageViewer();
      %viewer.setVisible(true);
      %imgs = javaArray('java.awt.image.BufferedImage',1)

      for nframe = frames

        % Load an image and convert it to the BufferedImage type of LOCI
        img = load_data(mymovie.(field)(k), nframe);

        if (crop_export)
          img = imalign(img, crop_size, mymovie.(lookup_field)(lookup_indx).centers(:, nframe+1), mymovie.(lookup_field)(lookup_indx).orientations(1, nframe+1));
        end

        img = im2java(img);
        img = AWTImageTools.makeBuffered(img, colorModel);

        %viewer.setImages(imgs);

        w.saveImage(img, 0, nframe == frames(end), nframe == frames(end));
      end

      try
        w.close();
      catch ME
        disp(ME.message);
        continue;
      end
    end
  end

  return;
end
