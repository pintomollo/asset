function [new_stack] = save_data(fname, data, varargin)

%% http://www.mathworks.com/matlabcentral/fileexchange/35684-save-and-load-data-as-multi-frame-tiff-format

  new_stack = [];
  if (nargin < 2)
    return;
  end

  types = {'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', 'uint64', 'int64'};
  split_rgb = false;
  metadata = '';
  type = '';
  for i=1:length(varargin)
    if (islogical(varargin{i}))
      split_rgb = varargin{i};
    elseif (ismember(varargin{i}, types))
      type = varargin{i};
    else
      metadata = varargin{i};
    end
  end

  if (~iscell(data))
    data = {data};
  end

  img = data{1};
  if (ischar(img))
    try
    img = imread(img);
    catch ME
      try
        img = movie2tiff(img);
      catch ME2
        error('ASSET:unknownFileFormat', ['The file ' img ' is of unknown format, please convert it manually into a stack of TIFF images.'])
      end
    end
  end

  is_partial=(nargout==1);

  is_sparse = issparse(img);

  [h,w,c] = size(img);
  split_rgb = (split_rgb && mod(c, 3)==0);

  if (isempty(type))
    type = class(img);
  end

  is_new = false;
  scaling_params = [];
  if (isstruct(fname))
    new_stack = fname;

    tiffobj = new_stack.tiffobj;
    tagstruct = new_stack.tagstruct;
    type = new_stack.type;
    files = new_stack.fname;
    scaling_params = new_stack.scaling;
  else
    if (type(1)=='u')
      byte = Tiff.SampleFormat.UInt;
    elseif (type(1)=='i')
      byte = Tiff.SampleFormat.Int;
    else
      warning('ASSET:saveData', ['Invalid pixel type ''' type '''. Using UINT8 instead.']);
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

    switch type
      case {'uint16', 'int16'}
          tagstruct.BitsPerSample = 16;
      case {'uint32', 'int32'}
          tagstruct.BitsPerSample = 32;
      case {'uint64', 'int64'}
          tagstruct.BitsPerSample = 64;
    end

    if (~isempty(metadata))
      tagstruct.ImageDescription = metadata;
    end

    is_new = true;
    if (split_rgb)
      if (iscell(fname))
        if (length(fname)<3)
          [file_path, filename, ext] = fileparts(fname{1});
          files = cell(3, 1);
          files{1} = fullfile(file_path, [filename '_R' ext]);
          files{2} = fullfile(file_path, [filename '_G' ext]);
          files{3} = fullfile(file_path, [filename '_B' ext]);
        else
          files = fname(1:3);
        end
      else
        [file_path, filename, ext] = fileparts(fname);
        files = cell(3, 1);
        files{1} = fullfile(file_path, [filename '_R' ext]);
        files{2} = fullfile(file_path, [filename '_G' ext]);
        files{3} = fullfile(file_path, [filename '_B' ext]);
      end
    else
      if (iscell(fname))
        files = fname(1);
      else
        files = {fname};
      end
    end
  end

  ncolors = 1 + 2*split_rgb;
  is_new = is_new(ones(1,ncolors));
  ndata = length(data);
  for i=1:ndata
    for n=1:c

      indx = mod(n, ncolors) + 1;

      if (~is_new(indx))
        tiffobj(indx).writeDirectory();
      else
        tiffobj(indx) = Tiff(files{indx}, 'w8');
        is_new(indx) = false;
      end
      tiffobj(indx).setTag(tagstruct);

      if (is_sparse)
        [cast_img, scaling_params] = scaled_cast(full(img(:,:,n)), scaling_params, type);
      else
        [cast_img, scaling_params] = scaled_cast(img(:,:,n), scaling_params, type);
      end
      tiffobj(indx).write(cast_img);
    end

    if (i~=ndata)
      img = data{i+1};
      if (ischar(img))
        img = imread(img);
      end

      is_sparse = issparse(img);

      [nh,nw,c] = size(img);

      if (nh ~= h || nw~=w || (split_rgb && mod(c, 3)~=0))
        for j=1:length(tiffobj)
          tiffobj(j).close();
        end
        error('Size of the provided images do not match !')
      end
    end
  end

  if (~is_partial)
    for j=1:length(tiffobj)
      tiffobj(j).close();
    end

    clearvars new_stack;
  else
    new_stack = struct('fname', files, ...
                       'tiffobj', tiffobj, ...
                       'tagstruct', tagstruct, ...
                       'scaling', scaling_params, ...
                       'type', type);
  end

  return;
end
