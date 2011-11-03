function [mymovie] = project3D(mymovie, proj_func, opts)

  if (nargin == 2)
    if (isstruct(proj_func))
      opts = proj_func;
      proj_func = @max;
    else
      opts = get_struct('ASSET');
    end
  end

  % Make sure the data have been copied to the current field
  if (~isfield(mymovie.data, 'centers') | isempty(mymovie.data.centers) | opts.recompute)
     mymovie = duplicate_segmentation(mymovie, 'data', opts);
  end

  % Make sure we have the metadata information
  if (~isfield(mymovie, 'metadata')|isempty(mymovie.metadata)|~isfield(mymovie.metadata, 'frame_index')|isempty(mymovie.metadata.frame_index))
    mymovie = parse_metadata(mymovie, opts);

    if (~isfield(mymovie, 'metadata')|isempty(mymovie.metadata)|~isfield(mymovie.metadata, 'frame_index')|isempty(mymovie.metadata.frame_index))
      warning('Sorry, cannot project a 4D recording without metadata !');

      return;
    end
  end

  channel_type = regexp(mymovie.data.file, '.+_(\w+)\.ome\.tiff?', 'tokens');

  if (isempty(channel_type))
    channel_index = 1;
  else
    channel_index = find(ismember(mymovie.metadata.channels, channel_type{1}), 1);

    if (isempty(channel_index))
      channel_index = 1;
    end
  end

  frames = mymovie.metadata.frame_index(channel_index, :);
  planes = mymovie.metadata.plane_index(channel_index, :);

  indexes = unique(frames);
  nframes = length(indexes);
  all_frames = [1:length(frames)];

  projection = get_struct('channel');
  ninputs = nargin(proj_func);
  fid = '';

  for nimg = indexes
    currents = (frames == nimg);
    stack = double(load_data(mymovie.data, all_frames(currents)));
    stack(:,:,planes(currents)) = stack;

    switch (ninputs)
      case 2
        projection_img = proj_func(stack, 3);
      case 3
        projection_img = proj_func(stack, [], 3);
      otherwise
        projection_img = proj_func(stack);
    end

    fid = save_data(fid, projection_img, nframes);
  end

  projection.fname = fid;
  mymovie.data.projection = projection;

  return;
end
