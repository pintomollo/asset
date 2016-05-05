function [new_stack] = save_stack(files, data, type, meta)

%% http://www.mathworks.com/matlabcentral/fileexchange/35684-save-and-load-data-as-multi-frame-tiff-format

  if (nargin > 1)
    if nargin==2, type = class(data); end

    if (nargin < 4), meta = ''; end

    if (iscell(files))
      files = files{1};
    end

    is_partial=(nargout==1);
    is_end = false;

    is_sparse = issparse(data);

    [h,w,c] = size(data);

    if (isstruct(files))
      new_stack = files;

      tiffobj = new_stack.tiffobj;
      tagstruct = new_stack.tagstruct;
      type = new_stack.type;

      ninit = 1;
      is_end = (~is_partial);
      is_partial = true;
    else

      if (type(1)=='u')
        byte = Tiff.SampleFormat.UInt;
      elseif (type(1)=='i')
        byte = Tiff.SampleFormat.Int;
      else
        warning(['Invalid pixel type ''' type '''. Using UINT8 instead.']);
        type = 'uint8';
        byte = Tiff.SampleFormat.UInt;
      end

      tagstruct = struct('ImageLength', h, ...
                         'ImageWidth', w, ...
                         'Photometric', Tiff.Photometric.MinIsBlack, ...
                         'BitsPerSample', 8, ...
                         'SamplesPerPixel', 1, ...
                         'RowsPerStrip', h, ...
                         'SampleFormat', byte, ...
                         'Compression', Tiff.Compression.None, ...
                         'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);

      switch class(data)
        case {'uint16', 'int16'}
            tagstruct.BitsPerSample = 16;
        case {'uint32', 'int32'}
            tagstruct.BitsPerSample = 32;
      end

      if (~isempty(meta))
        tagstruct.ImageDescription = meta;
      end

      if (~is_partial)
        msg = sprintf(' Writing planes of stack to disk : %3d', 1);
        fprintf(msg);
      end

      tiffobj = Tiff(files, 'w8');
      tiffobj.setTag(tagstruct);
      if (is_sparse)
        tiffobj.write(cast(full(data(:,:,1)), type));
      else
        tiffobj.write(cast(data(:,:,1), type));
      end

      new_stack = struct('fname', files, ...
                         'tiffobj', tiffobj, ...
                         'tagstruct', tagstruct, ...
                         'type', type);

      ninit = 2;
    end

    try
      for n = ninit:c % loop over images to resize images
        if (~is_partial)
          fprintf('\b\b\b%3d', n);
        end

        tiffobj.writeDirectory();
        tiffobj.setTag(tagstruct);

        if (is_sparse)
          tiffobj.write(cast(full(data(:,:,n)), type));
        else
          tiffobj.write(cast(data(:,:,n), type));
        end
      end
    catch exception
      tiffobj.close();

      throw(exception)
    end

    if (~is_partial)
      tiffobj.close();

      %fprintf('\b\b\b\bdone\n');
      fprintf([repmat('\b', 1, length(msg)) repmat(' ', 1, length(msg)) repmat('\b', 1, length(msg))]);
    elseif (is_end)
      tiffobj.close();
    end

  else

    if nargin<1, files = '*.tif'; end

    new_stack = '';
    dir_out = '_split'; % '../resized';

    [files, out_path] = get_filenames(files, dir_out);

    N = length(files);
    if (N == 0), disp('nada??'), return, end

    fprintf(' Reconstructing stack from planes :    1');

    filename = files{1};
    info = imfinfo(filename);
    im = imread(filename, 'Info', info);

    [h,w,c] = size(im);

    type = class(im);
    if (type(1)=='u')
      byte = Tiff.SampleFormat.UInt;
    else
      byte = Tiff.SampleFormat.Int;
    end
    tagstruct = struct('ImageLength', h, ...
                       'ImageWidth', w, ...
                       'Photometric', Tiff.Photometric.MinIsBlack, ...
                       'BitsPerSample', info.BitsPerSample(1), ...
                       'SamplesPerPixel', 1, ...
                       'RowsPerStrip', h, ...
                       'SampleFormat', byte, ...
                       'Compression', Tiff.Compression.None, ...
                       'PlanarConfiguration', Tiff.PlanarConfiguration.(info.PlanarConfiguration));

    new_stack = cell(c, 2);
    for i = 1:c
      new_stack{i,1} = fullfile(out_path, ['stack_' num2str(i) '.tif']);

      tiffobj = Tiff(new_stack{i,1}, 'w8');

      tiffobj.setTag(tagstruct);
      tiffobj.write(im(:,:,i));

      new_stack{i,2} = tiffobj;
    end

    try
      for n = 2:N % loop over images to resize images
        fprintf('\b\b\b%3d', n);

        filename = files{n};
        im = imread(filename);

        for i = 1:c
          tiffobj = new_stack{i,2};

          tiffobj.writeDirectory();
          tiffobj.setTag(tagstruct);
          tiffobj.write(im(:,:,i));
        end
      end
    catch exception
      for i = 1:c
        new_stack{i, 2}.close();
      end

      throw(exception)
    end

    for i = 1:c
      new_stack{i, 2}.close();
    end

    fprintf('\b\b\b\bdone\n');

    new_stack = new_stack(:,1);
  end

  return;
end
