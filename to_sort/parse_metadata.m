function [metadata, opts] = parse_metadata(fname, base_dir, opts)

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

    % Get the size of the problem, using the correct channel
    switch (opts.segmentation_type)

      % The default channel which we should always have
      case {'dic', 'all'}
        [nframes imgsize ] = size_data(mymovie.dic);

      % The marker segmentation is a bit more problematic as there might be no channel
      % for the eggshell
      case 'markers'
        if (isfield(mymovie, 'eggshell') & ~isempty(mymovie.eggshell))
          [nframes imgsize ] = size_data(mymovie.eggshell);
        else
          [nframes imgsize ] = size_data(mymovie.cortex);
        end

      case 'data'
        [nframes, imgsize] = size_data(mymovie.data);

      % Somehow the input parameters do not make sense
      otherwise
        error 'None of the expected field are present in ''mymovie''';
    end
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
    pixel_size = 0;
    binning = 0;
    magnification = 0;

    count = 1;
    state = 'init';
    block_level = 0;
    obj_level = 0;

    said_it = false;

    if (any(line == '{'))
      % uManager metadata file

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
                tmp_tokens = regexp(tokens{1}, 'FrameKey-(\d+)-\d+-\d+','tokens');
                if (~isempty(tmp_tokens) & str2double(tmp_tokens{1}{1}) >= nframes)
                  if (~said_it)
                    warning(['Metadata describes a plane (' tmp_tokens{1}{1}{1} ') that is out of the range of the actual recordings (' num2str(nframes-1) '), ignoring this infromation.']);
                    said_it = true;
                  end

                  state = 'ignore';
                  obj_level = block_level;
                else
                  state = 'read';
                  obj_level = block_level;
                end
              else
                state = 'ignore';
                obj_level = block_level;
              end
            end
          case 'ignore'
            if (~isempty(strfind(line, 'PixelSize')))
              tokens = regexp(line, '^\s*"(.+)": "?(.*?)"?,?$','tokens');
              pixel_size = str2double(tokens{1}{2});
            end
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

      if (~isempty(mymovie))
        if (all(plane(:) == 1))
          frame(:, end+1:nframes) = repmat([size(frame, 2)+1:nframes], size(frame, 1), 1);
          plane(:, end+1:nframes) = 1;
        else
          frame(:, end+1:nframes) = NaN;
          plane(:, end+1:nframes) = NaN;
        end
        pos(:, end+1:nframes) = NaN;
        time(:, end+1:nframes) = NaN;
        exposure(:, end+1:nframes) = NaN;
        group(:, end+1:nframes) = NaN;
      end

    elseif (any(line == '<'))
      % Spinning Disk metadata file

      % Loop until we have read the whole file
      while ischar(line)

        if (~isempty(strfind(line, '</')))
          if (sum(line == '<') == 1)
            block_level = block_level - 1;
          end
        elseif (isempty(strfind(line, '/>')))
          block_level = block_level + 1;
        end

        switch state
          case 'init'
            if (~any(line == '/'))
              tokens = regexp(line, '^\s*<(.+)>$','tokens');
              if (isempty(tokens))
                state = 'init';
              elseif (strncmp(tokens{1}, 'CameraSetting', 13)) | (strncmp(tokens{1}, 'SpatialCalibration', 18))
                state = 'read';
                obj_level = block_level;
              end
            end
          case 'read'
            if (block_level < obj_level)

              state = 'init';
            else
              tokens = regexp(line, '^\s*<(.+)>(.*?)</\1>$','tokens');
              if (~isempty(tokens))
                tokens = tokens{1};
              end
              
              if (numel(tokens) == 2)
                switch tokens{1}
                  case 'ObjectiveMagn'
                    magnification = str2double(tokens{2});
                  case 'XBasePixelSize'
                    pixel_size(1) = str2double(tokens{2});
                  case 'YBasePixelSize'
                    pixel_size(2) = str2double(tokens{2});
                  case 'Binning'
                    binning = str2double(tokens{2});
                  case 'ExposureTime'
                    exposure = str2double(tokens{2});
                  otherwise
                    %disp(['Parsing incomplete ...']);
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
    else
      warning('Unknown metadata file type, ignoring it.');

      return;
    end

    fclose(fid);

    metadata.z_position = pos;
    metadata.frame_index = frame;
    metadata.plane_index = plane;
    metadata.acquisition_time = time;
    metadata.exposure_time = exposure;
    metadata.channel_index = group;
    metadata.channels = groups;

    if (binning ~= 0 & magnification ~= 0 & pixel_size ~= 0)
      pixel_size = mean(pixel_size);

      opts.pixel_size = pixel_size;
      opts.magnification = magnification;
      opts.binning = binning;

      opts.ccd_pixel_size = pixel_size * magnification / binning;
    end
  end

  if (~isempty(mymovie))
    mymovie.metadata = metadata;
    metadata = mymovie;
  end

  return;
end
