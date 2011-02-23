function signal = quantify_signal(imgs, positions, opts)

  if (nargin == 2)
    if (isstruct(positions))
      opts = positions;
      positions = {};
    else
      opts = get_struct('ASSET');
    end
  elseif (nargin == 1)
    opts = get_struct('ASSET');
    positions = {};
  end

  channel_name = opts.quantification.channel;
  field_name = opts.quantification.field;

  if (isstruct(imgs))
    mymovie = imgs;
    imgs = [];

    if (opts.recompute | ~isfield(mymovie.(channel_name), 'eggshell') | isempty(mymovie.(channel_name).eggshell))
    
      mymovie = correct_dic_shift(mymovie, channel_name, opts.segmentation_parameters.correction, opts);
    end
    
    nframes = size_data(mymovie.(channel_name));
  else
    mymovie = [];
    nframes = size(imgs, 3);
  end

  window_shape = opts.quantification.window_shape;

  if (ischar(window_shape))
    window_size = opts.quantification.window_size / opts.pixel_size;

    switch window_shape
      case 'gaussian'
        window_shape = fspecial(window_shape, ceil(window_size), opts.quantification.window_params / opts.pixel_size);
        
      otherwise
        window_shape = fspecial(window_shape, ceil(window_size));
    end

    window_shape = window_shape / sum(window_shape(:));
  end

  for i=1:nframes
    if (isempty(mymovie))
      img = imgs(:,:,i);
      pos = positions{i};
    else
      img = load_data(mymovie.(channel_name), i);
      pos = mymovie.(channel_name).(field_name)(i).carth;
    end

    img = imfilter(img, window_shape);
    values = bilinear(img, pos(:,1), pos(:,2));
    
    if (isempty(imgs))
      mymovie.(channel_name).quantification(i).(field_name) = values;
    else
      signal{i} = values;
    end
  end

  if (isempty(imgs))
    signal = mymovie;
  end

  return;
end
