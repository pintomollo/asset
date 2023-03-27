function [metadata, opts] = parse_metadata(data, opts)
% PARSE_METADATA extracts relevant information from the metadata file.
%
%   [METADATA, OPTS] = PARSE_METADATA(XML, OPTS) tries to extract the acquisition
%   time, exposure time, z position and channel names of the recording and
%   stores them in METADATA. It also tries to extract the spatial resolution of
%   the image and stores it in OPTS.
%
%   [...] = PARSE_METADATA(UMANAGER, OPTS) converts a uManager text file into
%   XML before parsing it.
%
%   The XML formats currently accepted are :
%     - uManager text file     : 'www.w3.org/uManager'
%     - uManager metadata      : 'www.openmicroscopy.org/Schemas/OME/2012-06'
%     - PerkinElmer UltraVIEW  : 'tempuri.org/UltraviewSchm.xsd'
%     - Leica Application Suite: 'schemas.datacontract.org/2004/07/LeicaMicrosystems.DataEntities.V3_2'
%
% Wilson lab, University of Otago
% Simon Blanchoud
% 06.03.2015

  if (nargin < 2)
    opts = get_struct('options');
  end

  max_iter = 15;

  data = umanager2xml(data, max_iter);
  data = regexprep(data, '^\*.*\*', '');
  start = find(data(1:end-1)=='<' & data(2:end)=='?', 1, 'first');

  xml_data = parse_xml(data(start:end));

  xml_type = get_attribute(xml_data, 'xmlns');
  xml_type = regexp(xml_type, '^(http://)?(.*?)$', 'tokens');

  if (length(xml_type) > 0 && length(xml_type{1}) > 1)
    xml_type = xml_type{1}{2};

    switch xml_type
      case 'www.w3.org/uManager'
        summary_keys = {'Summary', 'Frames', 'Summary', 'Channels', 'Summary', 'Slices'};
        frame_keys = {'FrameKey*', '^Frame$', 'FrameKey*', '^Slice$', 'FrameKey*', '^Channel$', {'FrameKey*', 1/1000}, '^ElapsedTime-ms$', {'FrameKey*', 1/1000}, '^Exposure-ms$', {'FrameKey*', 1}, '^Z-um$'};
        infer_keys = {};
        resol_keys = {};
      case {'www.openmicroscopy.org/Schemas/OME/2012-06', 'www.openmicroscopy.org/Schemas/OME/2016-06'}
        summary_keys = {'Pixels', '^SizeT', 'Pixels', '^SizeC', 'Pixels', '^SizeZ'};
        frame_keys = {};
        infer_keys = {};
        resol_keys = {};
      case 'www.openmicroscopy.org/Schemas/OME/2015-01'
        summary_keys = {'Pixels', '^SizeT', 'Pixels', '^SizeC', 'Pixels', '^SizeZ'};
        frame_keys = {'Plane', 'TheT', 'Plane', 'TheZ', 'Plane', 'TheC', {'Plane', 1}, 'DeltaT', {'Plane', 1}, 'ExposureTime', {'Plane', 1}, 'PositionZ', 'UUID', 'FileName', 'Channel', 'Name'};
        infer_keys = {};
        resol_keys = {{'Pixels', 1}, {'PhysicalSizeX', 'PhysicalSizeY'}, '', '', {'Channel', 1}, 'SamplesPerPixel'};
      case 'tempuri.org/UltraviewSchm.xsd'
        summary_keys = {'Property', 'T', 'Property', 'C', 'Property', 'Z'};
        frame_keys = {};
        infer_keys = {'CameraSetting', 'ExposureTime','ChannelSetting', 'ChannelName',  {'AcquiredTime', 'yyyy-mm-ddTHH:MM:SS.FFF'}, {'StartTime', 'FinishTime'}, {'ZSetting', 1}, {'TopPosition', 'BottomPosition'}};
        resol_keys = {{'SpatialCalibration', 1}, {'XBasePixelSize', 'YBasePixelSize'}, {'SpatialCalibration', 1}, 'ObjectiveMagn', {'CameraSetting', 1}, 'Binning'};
      case 'schemas.datacontract.org/2004/07/LeicaMicrosystems.DataEntities.V3_2'
        summary_keys = {};
        frame_keys = {'', '', '', '', 'Camera', 'Name', {'LasImage', 'yyyy-mm-ddTHH:MM:SS.FFF'}, '^AcquiredDate$', {'Camera', '%f ms', 1/1000}, '^Exposure$', {'^Microscope_Focus_Position$', '%f mm', 1000}, ''};
        infer_keys = {};
        resol_keys = {{'LasImage', 1e6}, {'XMetresPerPixel', 'YMetresPerPixel'}, {'Microscope_Video_Magnification', 1}, '', '', ''};
      otherwise
        if (strncmp(xml_type, 'www.openmicroscopy.org/Schemas/OME/', 35))
          summary_keys = {'Pixels', '^SizeT', 'Pixels', '^SizeC', 'Pixels', '^SizeZ'};
          frame_keys = {};
          infer_keys = {};
          resol_keys = {};
          warning(['Attempting to parse an untested "' xml_type '" OME XML schema.']);
        else
          summary_keys = {};
          frame_keys = {};
          infer_keys = {};
          resol_keys = {};
          warning(['Unknown XML schema "' xml_type '", unable to parse it.']);
        end
    end

    metadata = parse_summary(xml_data, summary_keys);
    metadata = parse_frames(xml_data, frame_keys, metadata);
    metadata = infer_frames(xml_data, infer_keys, metadata);

    dt = diff(metadata.acquisition_time, [], 2);
    dt = median(dt(:));

    if (isfinite(dt) && dt > 0)
      opts.time_interval = dt;
    end

    metadata.raw_data = data;

    opts = infer_resolution(xml_data, resol_keys, opts);
  else
    metadata = get_struct('metadata');
    metadata = edit_options(metadata, 'Metadata missing, please provide some.');
  end

  return;
end

function data = convert_data(data, formats)

  if (isempty(formats) || isempty(data))
    return;
  end

  if (iscell(data))
    for i=1:numel(data)
      data{i} = convert_data(data{i}, formats);
    end
  else
    for i=1:numel(formats)
      format = formats{i};

      if (ischar(format))
        if (any(format == '%'))
          data = sscanf(data, format);
        else
          data = datevec(data, format);
        end
      elseif (isnumeric(format))
        if (ischar(data))
          data = str2double(data);
        end
        data = data * format;
      end
    end
  end

  return;
end

function value = get_attribute(node, attribute)

  value = '';
  if (isempty(node))
    return;
  end

  for i=1:length(node.Attributes)
    if (regexp(node.Attributes(i).Name, attribute))
      value = node.Attributes(i).Value;

      break;
    end
  end

  return;
end

function [channel_index, metadata] = get_channel(channel, metadata)

  is_empty = cellfun('isempty', metadata.channels);

  if (isempty(is_empty))
    channel_index = 1;
    metadata.channels{channel_index} = channel;
  else

    if (~all(is_empty))
      is_group = ismember(metadata.channels(~is_empty), channel);
    else
      is_group = ~is_empty;
    end

    if (any(is_group))
      channel_index = find(is_group, 1);
    else
      channel_index = find(is_empty, 1);

      metadata.channels{channel_index} = channel;
    end
  end

  return;
end

function child = get_child(node, children)

  child = '';

  for i=1:length(node.Children)
    if (regexp(node.Children(i).Name, children))
      if (isempty(child))
        child = node.Children(i);
      else
        child = [child node.Children(i)];
      end
    end
  end

  for i=1:length(node.Children)
    grand_child = get_child(node.Children(i), children);

    if (~isempty(grand_child))
      if (isempty(child))
        child = grand_child;
      else
        child = [child grand_child];
      end
    end
  end

  return;
end

function values = get_values(node, keys)

  nkeys = length(keys) / 2;

  values = cell(1, nkeys);
  has_any = false;

  for i=1:nkeys
    key = keys{2*i - 1};
    attr = keys{2*i};

    if (~isempty(key))
      formats = {};
      if (iscell(key))
        formats = key(2:end);
        key = key{1};
      end

      children = get_child(node, key);
      nchilds = length(children);

      if (nchilds > 0)
        good_child = false(1, nchilds);
        child_values = cell(nchilds, 1);

        for c = 1:nchilds
          child = children(c);

          if (isempty(attr))
            child_values{c} = child.Data;
          elseif (iscell(attr))
            nattrs = length(attr);
            tmp_val = cell(1, nattrs);

            for j=1:nattrs
              tmp_val{j} = get_attribute(child, attr{j});
            end
            child_values{c} = tmp_val;
          else
            child_values{c} = get_attribute(child, attr);
          end

          good_child(c) = good_child(c) || (~isempty(child_values{c}));
        end

        child_values = convert_data(child_values, formats);
      end

      child_values = child_values(good_child, :);
      values{i} = child_values;

      has_any = has_any || any(good_child);
    end
  end

  if (~has_any)
    values = {};
  end

  return;
end

function metadata = infer_frames(data, keys, metadata)

  if (isempty(keys))
    return;
  end

  [nchannels, nframes, nslices] = size(metadata.acquisition_time);

  values = get_values(data, keys);

  for i=1:nchannels
    if (~isempty(values{1}))
      metadata.exposure_time(i, :, :) = str2double(values{1}{i});
    end
    if (~isempty(values{2}))
      metadata.channels{i} = values{2}{i};
    end
  end

  if (~isempty(values{3}))
    if (numel(values{3}{1}{1}) > 1)
      start = values{3}{1}{1};
      ends = values{3}{1}{2};
    else
      start = datevec(values{3}{1}{1});
      ends = datevec(values{3}{1}{2});
    end
    step = etime(ends, start) / (nframes + 1);

    for i=1:nframes
      metadata.acquisition_time(:, i, :) = (i-1)*step;
    end
  end

  if (~isempty(values{4}))
    top = values{4}{1}{1};
    bottom = values{4}{1}{2};
    step = (top - bottom) / max(nslices-1, 1);

    for i=1:nslices
      metadata.z_position(:, :, i) = top + (i-1)*step;
    end
  end

  return;
end

function opts = infer_resolution(data, keys, opts)

  if (isempty(keys))
    return;
  end

  values = get_values(data, keys);

  if (isempty(values{1}))
    return;
  else
    nres = length(values{1}{1});
    pixel_size = NaN(1, nres);
    for i=1:nres
      pixel_size(i) = values{1}{1}{i};
    end
  end

  if (numel(pixel_size > 1) && all((pixel_size - pixel_size(1)) == 0))
    pixel_size = pixel_size(1);
  end

  if (isempty(values{2}))
    magnification = 1;
  else
    magnification = values{2}{1};
  end

  if (isempty(values{3}))
    binning = 1;
  else
    binning = values{3}{1};
  end

  if (binning > 0 && magnification > 0 && all(pixel_size > 0))

    opts.pixel_size = pixel_size;
    opts.magnification = magnification;
    opts.binning = binning;

    opts.ccd_pixel_size = pixel_size * magnification / binning;
  end

  return;
end

function metadata = parse_frames(data, keys, metadata)
  % Keys: 1. frames
  %       2. slice
  %       3. channel
  %       4. time
  %       5. exposure
  %       6. z-pos
  %       7. fname
  %       8. channel names

  if (isempty(keys))
    return;
  end

  tmp_channel = (round(rand(1,10)*57 + 65));
  tmp_channel(tmp_channel > 90 & tmp_channel < 97) = tmp_channel(tmp_channel > 90 & tmp_channel < 97) - 43;
  tmp_channel = char(tmp_channel);

  values = get_values(data, keys);
  nchild = max(cellfun('length', values));

  is_date = false;
  is_guess = false;

  count = 0;
  for i=1:nchild
    if (isempty(values{1}))
      frame = NaN;
    else
      frame = str2double(values{1}{i});
    end
    if (isempty(values{2}))
      slice = NaN;
    else
      slice = str2double(values{2}{i});
    end
    if (isempty(values{3}))
      channel = '';
    else
      channel = values{3}{i};
    end

    if (isnan(frame) && isnan(slice))
      is_guess = true;
      count = count + 1;
      frame = count;
      slice = 1;
    elseif (isnan(frame))
      frame = 1;
      slice = slice + 1;
    elseif (isnan(slice))
      frame = frame + 1;
      slice = 1;
    else
      frame = frame + 1;
      slice = slice + 1;
    end

    if (isempty(channel))
      channel = tmp_channel;
    end
    [channel, metadata] = get_channel(channel, metadata);

    if (isempty(values{4}))
      time = NaN;
    else
      time = values{4}{i};

      if (numel(time) > 1)
        time = datenum(time);
        is_date = true;
      end
    end
    if (isempty(values{5}))
      exposure = NaN;
    else
      exposure = values{5}{i};
    end
    if (isempty(values{6}))
      z_pos = NaN;
    else
      z_pos = values{6}{i};
    end
    if (isempty(values{7}))
      fname = '';
    else
      fname = values{7}{i};
    end

    metadata.acquisition_time(channel, frame, slice) = time;
    metadata.exposure_time(channel, frame, slice) = exposure;
    metadata.z_position(channel, frame, slice) = z_pos;
    metadata.files{channel, frame, slice} = fname;
  end

  valid_frames = ~any(any(isnan(metadata.acquisition_time),2), 1);
  if (any(valid_frames))

    metadata.acquisition_time = metadata.acquisition_time(:,:,valid_frames);
    metadata.exposure_time = metadata.exposure_time(:,:,valid_frames);
    metadata.z_position = metadata.z_position(:,:,valid_frames);
    metadata.files = metadata.files(:,:,valid_frames);

    if (is_guess)
      [junk, indexes] = sort(metadata.acquisition_time(1,1,:));

      metadata.acquisition_time = metadata.acquisition_time(:,:,indexes);
      metadata.exposure_time = metadata.exposure_time(:,:,indexes);
      metadata.z_position = metadata.z_position(:,:,indexes);
      metadata.files = metadata.files(:,:,indexes);
    end

    [junk, indexes] = sort(metadata.acquisition_time(:,1,1));

    metadata.channels = metadata.channels(indexes);
    metadata.acquisition_time = metadata.acquisition_time(indexes,:,:);
    metadata.exposure_time = metadata.exposure_time(indexes,:,:);
    metadata.z_position = metadata.z_position(indexes,:,:);
    metadata.files = metadata.files(indexes,:,:);

    metadata.acquisition_time = metadata.acquisition_time - metadata.acquisition_time(1);
    if (is_date)
      for i=1:numel(metadata.acquisition_time)
        metadata.acquisition_time(i) = etime(datevec(metadata.acquisition_time(i)), datevec(metadata.acquisition_time(1)));
      end
    end
  end

  %acq_times = metadata.acquisition_time(1,:,:);
  acq_times = metadata.acquisition_time;
  [nchannels, nframes, nslices] = size(acq_times);
  channels = repmat([1:nchannels].', 1, nframes, nslices);
  frames = repmat([1:nframes], nchannels, 1, nslices);
  slices = repmat(permute([1:nslices], [1 3 2]), nchannels, nframes, 1);

  [junk, indxs] = sort(acq_times(:));

  metadata.channel_index = channels(indxs).';
  metadata.frame_index = frames(indxs).';
  metadata.plane_index = slices(indxs).';

  if (length(values) > 7)
    if (length(metadata.channels) == length(values{8}))
      indexes = str2double(metadata.channels);

      if (min(indexes) == 0)
        indexes = indexes + 1;
      end

      metadata.channels = values{8}(indexes);
    end
  end

  return;
end

function metadata = parse_summary(xml_data, keys)

  metadata = get_struct('metadata');

  if (isempty(keys))
    return;
  end

  values = get_values(xml_data, keys);
  if (isempty(values))
    return;
  end

  nframes = str2double(values{1});
  nchannels = str2double(values{2});
  nslices = str2double(values{3});

  nframes = max(nframes, 1);
  nchannels = max(nchannels, 1);
  nslices = max(nslices, 1);

  data = NaN([nchannels, nframes, nslices]);

  metadata.channels = cell(nchannels, 1);

  metadata.acquisition_time = data;
  metadata.exposure_time = data;
  metadata.z_position = data;
  metadata.files = cell([nchannels, nslices, nframes]);

  return;
end

function data = umanager2xml(data, max_iter)

  if (any(data == '{'))
    rem_spaces = @remove_spaces;
    rem_eol = @remove_eol;

    data = regexprep(data, '("\S*? +\S*?":)', '${rem_spaces($1)}');
    data = regexprep(data, '(\[[^\]]*?\],)', '${rem_eol($1)}');

    data = regexprep(data, '"([^\n\r]*?)": ([^\n\r{]*?),?[\r\n]+', '<$1>$2</$1>\n');
    nprev = length(data);
    for i=1:max_iter
      data = regexprep(data, '[\n\r](\s+)"([^\n\r]*?)": {(.*?)[\n\r]\1}[,\r\n]+', '\n$1<$2>$3\n$1</$2>\n');
      curr = length(data);
      if (nprev == curr)
        break;
      else
        nprev = curr;
      end
    end
    data = regexprep(data, '{(.*)}', ['<?xml version="1.0" encoding="utf-8"?>\n'...
                                        '<TimeSeries xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '...
                                        'xmlns="http://www.w3.org/uManager">$1</TimeSeries>']);
    data = regexprep(data, '>"', '>');
    data = regexprep(data, '"<', '<');
  end

  return;
end

function str = remove_eol(str)

  str = regexprep(str, '\s*[\n\r]\s*', '');

  return;
end

function str = remove_spaces(str)

  str = strrep(str, ' ', '_');

  return;
end
