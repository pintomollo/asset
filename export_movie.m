function export_movie(mymovie, opts)

  if (nargin < 2)
    opts = get_struct('ASSET');
    %opts = struct('compression', '', ...
    %              'max_export', 1, ...
    %              'crop_export', true, ...
    %              'parse_export', 'normal');

  end

  %import loci.formats.ImageReader;
  %import loci.formats.out.OMETiffWriter;
  %import loci.formats.MetadataTools;

  %colorModel = AWTImageTools.makeColorModel(1, java.awt.image.DataBuffer.TYPE_BYTE);
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
          if (isfield(mymovie.(field)(k), 'centers') & ~isempty(mymovie.(field)(k).centers))
            lookup_field = field;
            lookup_indx = k;
          elseif (isfield(mymovie, 'markers') && isfield(mymovie.markers, 'centers'))
            lookup_field = 'markers';
            lookup_indx = 1;
          else
            crop_export = false;
            warning(['No segmentation data allows to crop the file ' mymovie.(field)(k).file]);
            continue;
          end
        end


        [tokens,junk]=regexp(mymovie.(field)(k).file,opts.file_regexpr,'tokens');
        name = tokens{1}{1};
        suffix = tokens{1}{2};
        ext = tokens{1}{3};

        if (~strcmp(ext,'.ome.tiff'))
          ext = '.ome.tiff';
        end

        if (crop_export)
          diff_name = 'CROP';
        else
          diff_name = 'ASSET';
        end

        if (~isempty(name))
          new_file = [name diff_name '-' suffix ext];
        else
          new_file = [suffix '-' diff_name ext];
        end

        [slash] = findstr(new_file,filesep);
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
              %try 
              %  r.close();
              %catch ME
              %  disp(ME.message);
              %  disp('Ignoring closing Error, Continuing');
              %end

              continue;
          end
        end
      end

      if (crop_export)

        if (length(opts.crop_size) == 1)
          crop_size = round(opts.crop_size * max(mymovie.(lookup_field)(lookup_indx).axes_length([2 1], :), [], 2));
        else
          crop_size = opts.crop_size;
          if (any(crop_size(:) < max(mymovie.(lookup_field)(lookup_indx).axes_length([2 1], :), [], 2)))
            warning('Crop size is smaller than some estimated embryos''s size')
          end
        end
      end
      %import loci.formats.gui.ImageViewer;

      if (do_init)
        [nframes, imgsize] = size_data(mymovie.(field)(k));

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

        frames = frames(1:max_frames);
        do_init = ~do_init;
      end

      % We convert an image to get the corresponding values for the metadata
      %if (isfield(mymovie.(field)(k), 'metadata'))
      %  omexmlMetaOrig = MetadataTools.createOMEXMLMetadata(mymovie.(field)(k).metadata);
%
%        new_xml = char(mymovie.(field)(k).metadata);
%        [starts, ends] = regexp(new_xml,'<Plane.*?>.*?</\Plane>');
%        new_xml = new_xml([1:ends(max_frames) (ends(end)+1):length(new_xml)]);
%
%        omexmlMeta = MetadataTools.createOMEXMLMetadata(java.lang.String(new_xml));
%        omexmlMeta.setPixelsPixelType(java.lang.String('uint8'), 0, 0);
%
%        has_metadata = true;
%      else
%        % could extend to export whole movie to a single file
%        channels = 1;
%
%        omexmlMeta = MetadataTools.createOMEXMLMetadata();
%        omexmlMeta.createRoot();
%
%        omexmlMeta.setPixelsSizeX(java.lang.Integer(imgsize(2)), 0, 0);
%        omexmlMeta.setPixelsSizeY(java.lang.Integer(imgsize(1)), 0, 0);
%        omexmlMeta.setPixelsSizeZ(java.lang.Integer(1), 0, 0);
%        omexmlMeta.setPixelsSizeC(java.lang.Integer(channels), 0, 0);
%        omexmlMeta.setPixelsPixelType(java.lang.String('uint8'), 0, 0);
%        omexmlMeta.setPixelsBigEndian(java.lang.Boolean(false), 0, 0);
%        omexmlMeta.setPixelsDimensionOrder(java.lang.String('XYCZT'), 0, 0);
%        omexmlMeta.setLogicalChannelSamplesPerPixel(java.lang.Integer(channels), 0, 0);
%
%        has_metadata = false;
%      end
    %omexmlMeta = loci.formats.MetadataTools.createOMEXMLMetadata();
    %r = loci.formats.ImageReader();
    %r.setMetadataStore(omexmlMeta);
    %r.setId(absolutepath(mymovie.(field)(k).file));

    %  if (crop_export)
    %    omexmlMeta.setPixelsSizeX(ome.xml.model.primitives.PositiveInteger(java.lang.Integer(crop_size(2))), 0);
    %    omexmlMeta.setPixelsSizeY(ome.xml.model.primitives.PositiveInteger(java.lang.Integer(crop_size(1))), 0);
%        omexmlMeta.setPixelsSizeX(java.lang.Integer(crop_size(2)), 0, 0);
%        omexmlMeta.setPixelsSizeY(java.lang.Integer(crop_size(1)), 0, 0);
    %  end

    %r.close();

%      omexmlMeta.setPixelsSizeT(java.lang.Integer(max_frames), 0, 0);
%      for i = 1:max_frames
%        nframe = frames(i);
%        if (has_metadata)
%          omexmlMeta.setPlaneTimingDeltaT(omexmlMetaOrig.getPlaneTimingDeltaT(0, 0, nframe), 0, 0, i - 1);
%        else
%          omexmlMeta.setPlaneTimingDeltaT(java.lang.Double(nframes * 10), 0, 0, i - 1);
%        end
%      end

      %if (isequal(new_file(1),'.'))
      %  new_file = fullfile(pwd, new_file(2:end));
      %end
      %new_file = absolutepath(new_file);

      %pType = char(omexmlMeta.getPixelsType(0).getValue());

      %vmin = double(intmin(pType));
      %vmax = double(intmax(pType));

      %keyboard

      %w = loci.formats.ImageWriter();
      %w = loci.formats.out.OMETiffWriter();
      %w = OMETiffWriter();
      %w = BufferedImageWriter(w);
      %w.setMetadataRetrieve(omexmlMeta);
      %w.setId(new_file);
      %w.setWriteSequentially(true);
      %w.setInterleaved(false);

      if(isempty(opts.compression))
        % Select the compression
        %compression = w.getCompressionTypes();
        compression = {'none', 'lzw', 'deflate', 'jpeg'};
        [sel,ok] = listdlg('PromptString','Select a Compression',...
                      'SelectionMode','single',...
                      'ListString',char(compression));
        if(ok==0)
          sel = 1;
        end

        opts.compression = compression(sel);
      end
      % Set the compression
      %w.setCompression(java.lang.String(opts.compression));

      %viewer = ImageViewer();
      %viewer.setVisible(true);
      %imgs = javaArray('java.awt.image.BufferedImage',1)

      nframes = length(frames);

      set(gca, 'position', [0 0 1 1], 'visible', 'off')

      for p = 1:nframes

        % Load an image and convert it to the BufferedImage type of LOCI
        %if (p == 1)
        %  [img, reader] = load_data(mymovie.(field)(k), frames(p));
        %else
          img = imnorm(double(load_data(mymovie.(field)(k), frames(p))));
        %end

        %keyboard

        if (crop_export)
          img = realign(img, crop_size, mymovie.(lookup_field)(lookup_indx).centers(:, frames(p)), mymovie.(lookup_field)(lookup_indx).orientations(1, frames(p)));
        end

        save_data(new_file, img);
        

        %{
        hold off;
        imshow(imnorm(img));
        hold on;

        switch field
          case 'dic'
            egg = realign(mymovie.(lookup_field)(lookup_indx).eggshell(frames(p)).carth, crop_size, mymovie.(lookup_field)(lookup_indx).centers(:, frames(p)), mymovie.(lookup_field)(lookup_indx).orientations(1, frames(p)));
            cortex = realign(mymovie.(lookup_field)(lookup_indx).cortex(frames(p)).carth, crop_size, mymovie.(lookup_field)(lookup_indx).centers(:, frames(p)), mymovie.(lookup_field)(lookup_indx).orientations(1, frames(p)));

            plot(egg(:,1), egg(:,2), 'Color', [93 255 107]/255, 'LineWidth', 1.6);
            plot(cortex(:,1), cortex(:,2), 'Color', [255 93 93]/255, 'LineWidth', 1.6);
          case 'cortex'
            cortex = insert_ruffles(mymovie.data.cortex(frames(p)).carth, mymovie.data.ruffles(frames(p)).paths);
            cortex = realign(cortex, crop_size, mymovie.(lookup_field)(lookup_indx).centers(:, frames(p)), mymovie.(lookup_field)(lookup_indx).orientations(1, frames(p)));
            ruffles = realign(mymovie.data.ruffles(frames(p)).carth, crop_size, mymovie.(lookup_field)(lookup_indx).centers(:, frames(p)), mymovie.(lookup_field)(lookup_indx).orientations(1, frames(p)));

            plot(cortex(:,1), cortex(:,2), 'Color', [255 93 93]/255, 'LineWidth', 1.6);
            scatter(ruffles(:,1), ruffles(:,2), 'MarkerEdgeColor', [93 101 255]/255, 'LineWidth', 1.5, 'SizeData', 16^2);
          case 'data'
            cortex = realign(mymovie.data.quantification(frames(p)).carth, crop_size, mymovie.(lookup_field)(lookup_indx).centers(:, frames(p)), mymovie.(lookup_field)(lookup_indx).orientations(1, frames(p)));

            plot(cortex(:,1), cortex(:,2), 'Color', [255 93 93]/255, 'LineWidth', 1.6);
        end

        print('-painters', '-dpng', ['PNG/' mymovie.experiment field '-' num2str(frames(p))  '.png']);
        %}

        %imwrite(img, new_file, 'TIFF', 'WriteMode', 'append');
        %new_file = save_data(new_file, img, nframes);

%        img = im2java(img);
%        img = AWTImageTools.makeBuffered(img, colorModel);

        %img = img * (vmax - vmin) + vmin;
        %img = cast(img, pType);

        %viewer.setImages(imgs);
        %img = img.';
        %img = 
        %img = typecast(img(:), 'uint8');
        %w.saveBytes(p - 1, img);
        %w.saveImage(img, 0, nframe == frames(end), nframe == frames(end));

      end

      %try
        %w.close();
        %reader.close();
      %catch ME
      %  disp(ME.message);
      %  continue;
      %end
    end
  end

  return;
end
