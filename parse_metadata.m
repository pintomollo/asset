function metadata = parse_metadata(fname, base_dir, opts)

  if (nargin < 3)
    if (isstruct(base_dir))
      opts = base_dir;
      base_dir = '.';
    else
      opts = get_struct('ASSET');
    end
  end

  metadata = get_struct('metadata');
  mymovie = [];
  z_pos = [];
  if (isstruct(fname))
    mymovie = fname;
    fname = mymovie.experiment;
  end

  if (~exist(fname, 'file'))
    if (exist([fname '.txt'], 'file') == 2)
      fname = [fname '.txt'];
    elseif (exist([base_dir filesep fname], 'file') == 2)
      fname = [base_dir filesep fname];
    elseif (exist([base_dir filesep fname '.txt'], 'file') == 2)
      fname = [base_dir filesep fname '.txt'];
    elseif (exist(['Movies' filesep fname], 'file') == 2)
      fname = ['Movies' filesep fname];
    elseif (exist(['Movies' filesep fname '.txt'], 'file') == 2)
      fname = ['Movies' filesep fname '.txt'];
    elseif (exist(['Movies' filesep 'metadata' filesep fname], 'file') == 2)
      fname = ['Movies' filesep 'metadata' filesep fname];
    elseif (exist(['Movies' filesep 'metadata' filesep fname '.txt'], 'file') == 2)
      fname = ['Movies' filesep 'metadata' filesep fname '.txt'];
    elseif (exist(['metadata' filesep fname], 'file') == 2)
      fname = ['metadata' filesep fname];
    elseif (exist(['metadata' filesep fname '.txt'], 'file') == 2)
      fname = ['metadata' filesep fname '.txt'];
    else
      fname = [];
    end
  end
  
  if (~isempty(fname))
    % Open the file, line by line
    fid = fopen(fname,'rt');

    % If it does not work, try the absolute path instead
    if (fid<0)
      fname = absolutepath(fname);
      fid = fopen(fname,'rt');
    end

    if (~isempty(mymovie) & ~isempty(mymovie.dic))
      z_pos = NaN(3, length(mymovie.dic.eggshell));
    end

    % Read the first line
    line = fgetl(fid);
    size_increment = 200;
    vector_size = size_increment;

    nans = NaN(1, size_increment);
    pos = nans;
    frame = nans;
    plane = nans;
    group = nans;
    exposure = nans;
    time = nans;
    groups = {};

    count = 1;
    state = 'init';
    block_level = 0;
    obj_level = 0;

    % Loop until we have read the whole file
    while ischar(line)
      if (any(line == '{'))
        block_level = block_level + 1;
      end
      if (any(line == '}'))
        block_level = block_level - 1;
      end

      switch state
        case 'init'
          if (any(line == '{'))
            tokens = regexp(line, '^\s*"(.+)": {$','tokens');
            if (isempty(tokens))
              state = 'init';
            elseif (strncmp(tokens{1}, 'FrameKey', 8))
              state = 'read';
              obj_level = block_level;
            else
              state = 'ignore';
              obj_level = block_level;
            end
          end
        case 'ignore'
          if (block_level < obj_level)
            state = 'init';
          end
        case 'read'
          if (block_level < obj_level)
            count = count + 1;
            if (count > vector_size)
              vector_size = vector_size + size_increment;
              pos = [pos nans];
              frame = [frame nans];
              plane = [plane nans];
              group = [group nans];
              exposure = [exposure nans];
              time = [time nans];
            end

            state = 'init';
          else
            tokens = regexp(line, '^\s*"(.+)": "?(.*?)"?,?$','tokens');
            if (~isempty(tokens))
              tokens = tokens{1};
            end
            
            if (numel(tokens) == 2)
              switch tokens{1}
                case 'Z-um'
                  pos(count) = str2double(tokens{2});
                case 'Frame'
                  frame(count) = str2double(tokens{2}) + 1;
                case 'Slice'
                  plane(count) = str2double(tokens{2}) + 1;
                case 'Channel'
                  is_group = ismember(groups, tokens{2});
                  if (isempty(is_group)|~any(is_group))
                    groups{end+1} = tokens{2};
                    is_group = [is_group true];
                  end
                  group(count) = find(is_group);
                case 'Exposure-ms'
                  exposure(count) = str2double(tokens{2});
                case 'ElapsedTime-ms'
                  time(count) = str2double(tokens{2});
                case {'Time', 'FileName'}
                  % Ignoring
                otherwise
                  disp(['Parsing incomplete ...']);
              end
            end
          end
        otherwise
          disp(['Unkown parsing state ("' state '"), skipping...']);
          state = 'init'; 
      end
    
      % Read a new line
      line = fgetl(fid);
    end

    fclose(fid);

    goods = (~isnan(frame) & ~isnan(plane));
    frame = frame(goods);
    plane = plane(goods);
    pos = pos(goods);
    time = time(goods);
    exposure = exposure(goods);
    group = group(goods);

    frames = unique(frame);
    planes = unique(plane);
    channels = unique(group);
    sizes = [length(channels), length(frames), length(planes)];

    order = sub2ind(sizes(2:3), frame, plane);
    sizes = [sizes(1) prod(sizes(2:3))];
    order = sub2ind(sizes, group, order);

    tmp = NaN(sizes);
    tmp(order) = pos;
    pos = tmp;
    tmp(order) = frame;
    frame = tmp;
    tmp(order) = plane;
    plane = tmp;
    tmp(order) = time;
    time = tmp;
    tmp(order) = exposure;
    exposure = tmp;
    tmp(order) = group;
    group = tmp;

    [iindx, jindx] = find(isnan(pos));

    for i = iindx
      for j = jindx
        if (j > 1)
          pos(i, j) = pos(i, j-1);
          frame(i, j) = frame(i, j-1) + 1;
          plane(i, j) = plane(i, j-1);
          time(i, j) = time(i, j-1);
          exposure(i, j) = exposure(i, j-1);
          group(i, j) = group(i, j-1);
        end
      end
    end

    metadata.z_position = pos;
    metadata.frame_index = frame;
    metadata.plane_index = plane;
    metadata.acquisition_time = time;
    metadata.exposure_time = exposure;
    metadata.channel_index = group;
    metadata.channels = groups;
  end

  if (~isempty(mymovie))
    mymovie.metadata = metadata;
    metadata = mymovie;
  end

  return;
end
