function [myrecording, opts, is_updated] = inspect_recording(fname, opts)
% INSPECT_RECORDING displays a pop-up window for the user to manually identify the
% type of data contained in the different channels of a movie recording.
%
%   [MYRECORDING] = INSPECT_RECORDING(CHANNELS) displays the window using the data
%   contained in CHANNELS, updates it accordingly to the user's choice and returns
%   the adequate structure for later analysis MYRECORDING. CHANNELS can either
%   be a string, a cell list of strings or a 'channel' structure (see get_struct.m).
%   MYRECORDING is a structure as defined by get_struct('myrecording').
%
%   [...] = INSPECT_RECORDING() prompts the user to select a recording and converts
%   it before opening the GUI.
%
%   [MYRECORDING, OPTS] = INSPECT_RECORDING(...) returns in addition the parameters
%   required to filter the various channels as chosen by the user.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 14.05.2014

  % Argument checking, need to know if we ask for a recording or not.
  if (nargin == 0 || isempty(fname) || ...
     (isstruct(fname) && isfield(fname, 'channels') && isempty(fname.channels)))
    fname = convert_movie();
  end

  % We did not get anything to handle...
  if isempty(fname)
    myrecording = [];
    opts = [];
    is_updated = false;

    return;
  end

  % The structure containing the parameters for the different filters available to
  % the user
  if (nargin < 2 || isempty(opts))
    opts = get_struct('options');
  end

  % Store the original options
  orig_opts = opts;

  % Was it a tracking file ?
  was_tracking = false;

  % Create the channels structure if it was not provided.
  if (isstruct(fname))
    if (isfield(fname, 'experiment'))
      myrecording = fname;
      channels = myrecording.channels;
      was_tracking = true;

      % Retrieve the original file
      for i = 1:length(channels)
        if (~isempty(channels(i).file))
          channels(i).fname = channels(i).file;
        end
      end
    else
      channels = fname;
    end
  else
    % Put everything in a cell list
    if (ischar(fname))
      fname = {fname};
    end

    % Parse the list and copy its content into the channels structure
    nchannels = length(fname);
    channels = get_struct('channel', [nchannels 1]);
    for i=1:nchannels
      channels(i).fname = fname{i};
    end
  end

  % Dragzoom help message
  imghelp = regexp(help('dragzoom'), ...
             '([ ]+Normal mode:.*\S)\s+Mouse actions in 3D','tokens');
  imghelp = ['DRAGZOOM interactions (help dragzoom):\n\n', imghelp{1}{1}];

  % Create the GUI
  [hFig, handles] = create_figure();

  % Allocate various variables. This allows them to be "persistent" between
  % different calls to the callback functions.
  img = [];
  orig_img = [];
  img_next = [];
  is_updated = true;

  % And handle the colormaps as well
  colors = get_struct('colors');
  color_index = 1;

  % Display the figure
  set(hFig,'Visible', 'on');
  % Update its content
  update_display;
  % And wait until the user is done
  uiwait(hFig);

  % Now that the data are correct, create the whole structure
  if (~was_tracking || is_updated)
    myrecording = get_struct('myrecording');
  end

  if (is_updated)
    % Copy the channels
    myrecording.channels = channels;
    % And get the experiment name
    myrecording.experiment = get(handles.experiment, 'String');
  end

  % Delete the whole figure
  delete(hFig);
  drawnow;

  return;

  function update_display(recompute)
  % The main figure of the GUI, the one responsible for the proper display
  % of its content.

    % By default we recompute everything
    if (nargin == 0)
      recompute = true;
    end

    % Get the indexes of the current frame and channel
    indx = handles.current;
    nimg = handles.frame;

    % Stop if no data at all
    if (indx < 1)
      return;
    end

    % If we have changed channel, we need to update the display of the buttons
    if (indx ~= handles.prev_channel)

      % Get the colormap for the displayed channel
      color_index = channels(indx).color;

      % Set the name of the current panel
      set(handles.uipanel,'Title', ['Channel ' num2str(indx)]);

      % The filters
      set(handles.detrend,'Value', channels(indx).detrend);
      set(handles.hot_pixels,'Value', channels(indx).hot_pixels);
      set(handles.normalize,'Value', channels(indx).normalize);
      set(handles.cosmics,'Value', channels(indx).cosmics);

      % The type and compression
      set(handles.channel_type, 'Value',  channels(indx).type);
      set(handles.compress, 'Value',  channels(indx).compression);

      % And setup the indexes correctly
      handles.prev_channel = indx;
      handles.prev_frame = -1;
    end

    % Here we recompute all the filtering of the frame
    if (recompute)
      % Because it takes long, display it and block the GUI
      set(hFig, 'Name', 'Channel Identification (Filtering...)');
      set(handles.text, 'String', ['Frame #' num2str(nimg)]);
      set(handles.all_buttons, 'Enable', 'off');
      drawnow;
      refresh(hFig);

      % Try to avoid reloading frames as much as possible
      if (handles.prev_frame == nimg-1)
        orig_img = img_next;
        img_next = double(load_data(channels(indx).fname, nimg+1));
      elseif (handles.prev_frame == nimg+1)
        img_next = orig_img;
        orig_img = double(load_data(channels(indx).fname, nimg));
      elseif (handles.prev_frame ~= nimg)
        orig_img = double(load_data(channels(indx).fname, nimg));
        img_next = double(load_data(channels(indx).fname, nimg+1));
      end

      % Update the index
      handles.prev_frame = nimg;

      % Copy to the working variable
      img = orig_img;

      % Detrend the image ?
      if (channels(indx).detrend)
        img = imdetrend(img, opts.filtering.detrend_meshpoints);
      end

      % Remove cosmic rays ?
      if (channels(indx).cosmics)
        img = imcosmics(img, opts.filtering.cosmic_rays_window_size, opts.filtering.cosmic_rays_threshold);
      end

      % Remove hot pixels ?
      if (channels(indx).hot_pixels)
        img = imhotpixels(img, opts.filtering.hot_pixels_threshold);
      end

      % Normalize the image ?
      if (channels(indx).normalize)
        img = imnorm(img);
      end
    end

    % Determine which image to display in the left panel
    switch handles.display(1)

      % The raw image
      case 2
        img1 = orig_img;

      % The difference between filtered and raw
      case 3
        if (channels(indx).normalize)
          img1 = (imnorm(orig_img) - img);
        else
          img1 = orig_img - img;
        end

      % The filtered image
      otherwise
        img1 = img;
    end

    % Determine which image to display in the right panel
    switch handles.display(2)

      % The next frame
      case 2
        img2 = img_next;

      % The difference between current and next frame
      case 3
        if (isempty(img_next))
          img2 = [];
        else
          img2 = (orig_img - img_next);
        end

      % The raw image
      otherwise
        img2 = orig_img;
    end

    % If we have already created the axes and the images, we can simply change their
    % content (i.e. CData)
    if (numel(handles.img) > 1 & all(ishandle(handles.img)))
      set(handles.img(1),'CData', img1);
      set(handles.img(2),'CData', img2);
    else

      % Otherwise, we create the two images in their respective axes
      handles.img = image(img1,'Parent', handles.axes(1),...
                        'CDataMapping', 'scaled',...
                        'Tag', 'image');
      handles.img(2) = image(img2,'Parent', handles.axes(2), ...
                        'CDataMapping', 'scaled',...
                        'Tag', 'image');

      % Hide the axes and prevent a distortion of the image due to stretching
      set(handles.axes,'Visible', 'off',  ...
                 'DataAspectRatio',  [1 1 1]);

      % Drag and Zoom library from Evgeny Pr aka iroln
      dragzoom(handles.axes, 'on')
    end

    % And set the colormap
    if (ischar(colors.colormaps{color_index}))
      colormap(hFig, brewermap(128, colors.colormaps{color_index}));
    else
      colormap(hFig, colors.colormaps{color_index}(128));
    end

    if (recompute)
      % Release the image
      set(hFig, 'Name', 'Channel Identification');
      set(handles.all_buttons, 'Enable', 'on');
    end

    return
  end

  function remove_channel_Callback(hObject, eventdata)
  % This function removes ones channel from the list

    % Get the current index for later
    indx = handles.current;

    % And ask for confirmation
    answer = questdlg(['Are you sure you want to remove channel ' num2str(indx) ' ?' ...
                '(No data will be deleted from the disk)'], 'Removing a channel ?');
    ok = strcmp(answer, 'Yes');

    % If it's ok, let's go
    if (ok)

      % Remove the current index and select the first one
      channels(indx) = [];
      handles.current = 1;

      % If it was the only one, we need to handle this
      if (isempty(channels))

        % Set up the indexes as empty
        handles.current = 0;
        handles.prev_channel = -1;

        % As this is long, block the GUI
        set(handles.all_buttons, 'Enable', 'off');
        set(handles.img, 'CData', []);
        set(handles.list, 'String', '');

        % And call the movie conversion function to load a new one
        new_channel = convert_movie();

        % Did we get something back ?
        if (~isempty(new_channel))

          % Get the number of frames in each of them
          nframes = size_data(new_channel);

          % Handle single frame
          set(handles.slider, 'Max', max(nframes,1.1), 'Value', 1);

          % We provide basic default values for all fields
          channels = get_struct('channel');
          channels.fname = new_channel;
          channels.type = 1;
          channels.compression = 1;

          % We update the list of available channels
          set(handles.list, 'String', 'Channel 1', 'Value', 1);
        end

      % Otherwise, delete the last traces of the deleted channel
      else
        handles.prev_channel = -1;
        tmp_list = get(handles.list, 'String');
        tmp_list(indx,:) = [];
        set(handles.list, 'String', tmp_list, 'Value', 1);
      end

      % Release the GUI and update
      set(handles.all_buttons, 'Enable', 'on');
      update_display(true);
    end

    return;
  end

  function add_channel_Callback(hObject, eventdata)
  % This function adds a new channel to the current recording

    % Remember how many there were
    nchannels = length(channels);

    % As this is long, block the GUI
    set(handles.all_buttons, 'Enable', 'off');

    % And call the movie conversion function to get a new channel
    new_channel = convert_movie();

    % Did we get anything ?
    if (~isempty(new_channel))

      % Get the number of frames in each of them
      nframes = size_data(new_channel);
      curr_nframes = round(get(handles.slider, 'Max'));

      % If they are similar, we can add it to the current structure
      if (nchannels == 0 || nframes == curr_nframes)

        % We provide basic default values for all fields
        channels(end+1) = get_struct('channel');
        channels(end).fname = new_channel;
        channels(end).type = 1;
        channels(end).compression = 1;

        % We update the list of available channels
        tmp_list = get(handles.list, 'String');
        if (nchannels == 0)
          liststring = ['Channel ' num2str(length(channels))];
        else
          liststring = [tmp_list '|Channel ' num2str(length(channels))];
        end
        set(handles.list, 'String', liststring);

      % Otherwise, there is a problem !
      else
        errordlg(['Error: the selected channel does not have the same number of ' ...
                ' frames (' num2str(nframes) ') than the currently loaded ones (' ...
                num2str(curr_nframes) '), ignoring it.'],'Error: Adding a channel','modal');
      end
    end

    % Release the GUI
    set(handles.all_buttons, 'Enable', 'on');

    % If we just add a new first channel, we need to set it up properly and update
    if (nchannels == 0 && length(channels) > 0)
      handles.current = 1;
      set(handles.list, 'Value', 1);
      update_display(true);
    end

    return
  end

  function options_Callback(hObject, eventdata)
  % This function is responsible for handling the buttons responsible for the
  % option structure

    % Block the GUI
    set(handles.all_buttons, 'Enable', 'off');
    drawnow;
    refresh(hFig);

    % And get the type of button which called the callback (from its tag)
    type = get(hObject, 'tag');

    % By default, recompute
    recompute = true;

    % Handle all three buttons differently
    switch type

      % Call the editing function
      case 'edit'
        [opts.filtering, recompute] = edit_options(opts.filtering);

      % Call the loading function
      case 'load'
        opts = load_parameters(opts);

      % Call the saving function
      case 'save'
        save_parameters(opts);
        recompute = false;

      % Save a snapshot
      case 'snapshot'

        % Fancy output
        disp('[Select a SVG filename]');

        % Prompting the user for the filename
        [fname, dirpath] = uiputfile({'*.svg', 'SVG vectorized image'}, ['Select a filename for your snapshot'], 'export/snapshot.svg');

        % Not cancelled
        if (ischar(fname))

          % This might take a while
          curr_name = get(hFig, 'Name');
          set(hFig, 'Name', [curr_name ' (Saving snapshot...)']);

          % Get the full name and save the snapshot !
          fname = fullfile(dirpath, fname);
          plot2svg(fname, hFig);

          % And release !
          set(hFig, 'Name', curr_name);
        end

        recompute = false;
    end

    % Release the GUI and recompute the display
    set(handles.all_buttons, 'Enable', 'on');
    update_display(recompute);

    return
  end

  function gui_Callback(hObject, eventdata)
  % This function handles the callback of most buttons in the GUI !

    % By default we recompute the display
    recompute = true;

    % Get the channel index
    indx = handles.current;

    % If no more data, do nothing
    if (indx < 1)
      return;
    end

    % And get the type of button which called the callback (from its tag)
    type = get(hObject, 'tag');
    switch type

      % Each checkbox is responsible for its respective boolean fields
      case {'detrend', 'cosmics', 'hot_pixels', 'normalize'}
        channels(indx).(type) = logical(get(hObject, 'Value'));

      % The slider varies the frame index
      case 'slider'
        handles.frame = round(get(hObject, 'Value'));

      % The radio buttons have the index of their respective choice encoded
      % in their tag (e.g. radioXY). However, because all the images are stored
      % we do not need to recompute anything !
      case 'radio'
        tmp_tag = get(eventdata.NewValue, 'tag');
        handles.display(str2double(tmp_tag(end-1))) = [str2double(tmp_tag(end))];
        recompute = false;

      % A change in the channel index
      case 'channels'
        handles.current = get(hObject, 'Value');

      % A different selection in one of the drop-down lists
      case {'type', 'compression'}
        channels(indx).(type) = get(hObject, 'Value');
        recompute = false;

      % Call the color gui
      case 'color'
        [tmp_index, recompute] = gui_colors(color_index);
        if (recompute)
          color_index = tmp_index;
          channels(indx).color = color_index;
        end

      % Otherwise, do nothing. This is used to cancel the deletion requests
      otherwise
        return;
    end

    % Update the display accordingly
    update_display(recompute);

    return
  end

  function cancel_CloseRequestFcn(hObject, eventdata)
  % This function stops the current processing after confirmation

    % Just double check that the user want to quit
    answer = questdlg('Do you really want to discard all your changes ?');
    ok = strcmp(answer,'Yes');

    % If everything is OK, release the GUI and quit
    if (ok)
      is_updated = false;
      opts = orig_opts;
      uiresume(hFig);
    end

    return
  end

  function channel_CloseRequestFcn(hObject, eventdata)
  % This function converts the various indexes back into strings to prepare
  % the channels structure for its standard form before releasing the GUI
  % to exit it

    % Create a copy of channels in case we need to cancel
    tmp_channels = channels;

    % We need ot loop over the channels
    nchannels = length(channels);

    % Get the available types of channels and compressions
    contents = get(handles.channel_type,'String');
    ntypes = length(contents);
    compressions = get(handles.compress,'String');

    % We want to check if some channels have identical features
    detrend = logical(zeros(nchannels,1));
    types = logical(zeros(nchannels,ntypes));

    % Convert the indexes into strings and build a summary of the filters
    for i=1:nchannels
      detrend(i) = channels(i).detrend;
      types(i,channels(i).type) = true;
      tmp_channels(i).type = contents{channels(i).type};
      tmp_channels(i).compression = compressions{channels(i).compression};
    end

    % Some checks to make sure the user is aware of some potential issues
    ok = true;
    if (ok & any(detrend, 1))

      % This is a highly non-linear filtering which prevents proper
      % signal comparison between recordings
      answer = questdlg({'Some channels will be detrended, continue ?','', ...
                         '(This is a quite slow and non-linear process..)'});
      ok = strcmp(answer,'Yes');
    end

    % If everything is OK, release the GUI and quit
    if (ok)
      channels = tmp_channels;
      uiresume(hFig);
    end

    return
  end

  function [hFig, handles] = create_figure
  % This function actually creates the GUI, placing all the elements
  % and linking the callbacks.

    % The number of channels provided
    nchannels = length(channels);

    % Initialize the possible types and compressions
    typestring = {'brightfield'; 'dic'; 'eggshell'; 'cortex'; 'fluorescence'};
    typecompress = {'none', 'lzw', 'deflate', 'jpeg'};

    % Initialize the structure used for the interface
    liststring = '';
    for i = 1:nchannels

      % Build the displayed list
      liststring = [liststring 'Channel ' num2str(i)];
      if (i < nchannels)
        liststring = [liststring '|'];
      end

      % Set the currently selected type of data
      for j = 1:length(typestring)
        if (strcmp(channels(i).type, typestring{j}))
          channels(i).type = j;

          break;
        end
      end

      % If the type did not exist, use the first one
      if (ischar(channels(i).type))
        channels(i).type = 1
      end

      % Set the compression type
      for j = 1:length(typecompress)
        if (strcmpi(channels(i).compression, typecompress{j}))
          channels(i).compression = j;

          break;
        end
      end

      % If none was found, choose the first one
      if (ischar(channels(i).compression))
        channels(i).compress = 1;
      end
    end

    % Get the number of frames
    nframes = size_data(channels(1).fname);

    % Create a name for the experiment based on the filename
    exp_name = channels(1).fname;
    for i = 2:nchannels
      exp_name = common_substring(exp_name, channels(i).fname);
    end
    [junk, exp_name, junk] = fileparts(exp_name);
    [junk, exp_name, junk] = fileparts(exp_name);
    exp_name = regexprep(exp_name, ' ', '');

    % Create my own grayscale map for the image display
    mygray = [0:255]' / 255;
    mygray = [mygray mygray mygray];

    % We build a list of all buttons to easily block and release them
    enabled = [];

    % The main figure, cannot be rescaled, closed nor deleted
    hFig = figure('PaperUnits', 'centimeters',  ...
                  'CloseRequestFcn', @cancel_CloseRequestFcn, ...
                  'Color',  [0.7 0.7 0.7], ...
                  'Colormap', mygray, ...
                  'MenuBar', 'none',  ...
                  'Name', 'Channel Identification',  ...
                  'NumberTitle', 'off',  ...
                  'Units', 'normalized', ...
                  'Position', [0 0 1 1], ...
                  'DeleteFcn', @gui_Callback, ...
                  'HandleVisibility', 'callback',  ...
                  'Tag', 'channel_fig',  ...
                  'UserData', [], ...
                  'Visible', 'off');

    %%%%%% Now the buttons around the main panel

    % The list of channels
    hChannel = uicontrol('Parent', hFig, ...
                         'Units', 'normalized',  ...
                         'Callback', @gui_Callback, ...
                         'Position', [0.01 0.11 0.1 0.79], ...
                         'String', liststring, ...
                         'Style', 'listbox',  ...
                         'Value', 1, ...
                         'Tag', 'channels');
    enabled = [enabled hChannel];

    % The OK button
    hOK = uicontrol('Parent', hFig, ...
                    'Units', 'normalized',  ...
                    'Callback', @channel_CloseRequestFcn, ...
                    'Position', [0.70 0.02 0.18 0.05], ...
                    'String', 'OK',  ...
                    'Tag', 'pushbutton11');
    enabled = [enabled hOK];

    % The Cancel button
    hCancel = uicontrol('Parent', hFig, ...
                    'Units', 'normalized',  ...
                    'Callback', @cancel_CloseRequestFcn, ...
                    'Position', [0.90 0.02 0.08 0.05], ...
                    'String', 'Cancel',  ...
                    'Tag', 'pushbutton12');
    enabled = [enabled hCancel];

    % The Add and Remove buttons
    hAdd = uicontrol('Parent', hFig, ...
                    'Units', 'normalized',  ...
                    'Callback', @add_channel_Callback, ...
                    'Position', [0.01 0.055 0.1 0.04], ...
                    'String', 'Add channel',  ...
                    'Tag', 'pushbutton13');
    enabled = [enabled hAdd];

    hRemove = uicontrol('Parent', hFig, ...
                    'Units', 'normalized',  ...
                    'Callback', @remove_channel_Callback, ...
                    'Position', [0.01 0.01 0.1 0.04], ...
                    'String', 'Remove channel',  ...
                    'Tag', 'pushbutton14');
    enabled = [enabled hRemove];

    % The experiment name and its labels
    hText = uicontrol('Parent', hFig, ...
                      'Units', 'normalized',  ...
                      'Position', [0.2 0.93 0.09 0.025], ...
                      'String', 'Experiment name:',  ...
                      'TooltipString', sprintf(imghelp), ...
                      'BackgroundColor', get(hFig, 'Color'), ...
                      'FontSize', 12, ...
                      'Style', 'text',  ...
                      'Tag', 'text1');

    hName = uicontrol('Parent', hFig, ...
                      'Units', 'normalized',  ...
                      'Position', [0.3 0.93 0.5 0.05], ...
                      'String', exp_name,  ...
                      'FontSize', 12, ...
                      'Style', 'edit',  ...
                      'Tag', 'experiment');
    enabled = [enabled hName];

    % The slider and its label
    hIndex = uicontrol('Parent', hFig, ...
                      'Units', 'normalized',  ...
                      'Position', [0.2 0.03 0.09 0.025], ...
                      'String', 'Frame #1',  ...
                      'BackgroundColor', get(hFig, 'Color'), ...
                      'FontSize', 12, ...
                      'Style', 'text',  ...
                      'Tag', 'text2');

    hFrame = uicontrol('Parent', hFig, ...
                    'Units', 'normalized',  ...
                    'Callback', @gui_Callback, ...
                    'Position', [0.3 0.03 0.35 0.025], ...
                    'Value', 1, ...
                    'SliderStep', [1 10]/nframes, ...
                    'Max', max(nframes, 1.1), ...
                    'Min', 1, ...
                    'Style', 'slider', ...
                    'Tag', 'slider');
    enabled = [enabled hFrame];

    %%%%%%% Now the main panel

    % The panel itsel
    hPanel = uipanel('Parent', hFig, ...
                     'Title', 'Channel 1',  ...
                     'Tag', 'uipanel',  ...
                     'Clipping', 'on',  ...
                     'Position', [0.12 0.11 0.87 0.8]);

    % The two axes
    hAxes = axes('Parent', hPanel, ...
                 'Position', [0 0.1 0.43 0.9], ...
                 'DataAspectRatio', [1 1 1], ...
                 'Visible', 'off',  ...
                 'Tag', 'axes');

    hAxesNext = axes('Parent', hPanel, ...
                 'Position', [0.44 0.1 0.43 0.9], ...
                 'DataAspectRatio', [1 1 1], ...
                 'Visible', 'off',  ...
                 'Tag', 'axes');

    % The two radio button groups that handle which image to display
    % For the choices to be mutually exclusive, one has to put them inside
    % such uibuttongroup.
    hRadio = uibuttongroup('Parent', hPanel, ...
                         'Units', 'normalized',  ...
                         'SelectionChangeFcn', @gui_Callback, ...
                         'Position', [0 0.05 0.43 0.05], ...
                         'tag', 'radio');

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.1 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Filtered image', ...
                         'Tag', 'radio11');
    enabled = [enabled hControl];

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.4 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Raw image', ...
                         'Tag', 'radio12');
    enabled = [enabled hControl];

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.7 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Difference', ...
                         'Tag', 'radio13');
    enabled = [enabled hControl];

    % The second group
    hRadio = uibuttongroup('Parent', hPanel, ...
                         'Units', 'normalized',  ...
                         'SelectionChangeFcn', @gui_Callback, ...
                         'Position', [0.44 0.05 0.43 0.05], ...
                         'tag', 'radio');

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.1 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Current frame', ...
                         'Tag', 'radio21');
    enabled = [enabled hControl];

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.4 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Next frame', ...
                         'Tag', 'radio22');
    enabled = [enabled hControl];

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.7 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Difference', ...
                         'Tag', 'radio23');
    enabled = [enabled hControl];

    % The type, color and compression of the channel, along with its labels
    hText = uicontrol('Parent', hPanel, ...
                      'Units', 'normalized',  ...
                      'Position', [0.9 0.925 0.05 0.05], ...
                      'String', 'Channel:',  ...
                      'FontSize', 12, ...
                      'FontWeight', 'bold', ...
                      'Style', 'text',  ...
                      'Tag', 'text16');

    hText = uicontrol('Parent', hPanel, ...
                      'Units', 'normalized',  ...
                      'Position', [0.9 0.875 0.05 0.05], ...
                      'String', 'Type',  ...
                      'Style', 'text',  ...
                      'Tag', 'text17');

    hType = uicontrol('Parent', hPanel, ...
                      'Units', 'normalized',  ...
                      'Callback', @gui_Callback, ...
                      'Position', [0.875 0.845 0.1 0.05], ...
                      'String', typestring, ...
                      'Style', 'popupmenu',  ...
                      'Value', 1, ...
                      'Tag', 'type');
    enabled = [enabled hType];

    % The buttons which allows to change the colormap
    hColor = uicontrol('Parent', hPanel, ...
                       'Units', 'normalized',  ...
                       'Callback', @gui_Callback, ...
                       'Position', [0.89 0.77 0.08 0.04], ...
                       'Style', 'pushbutton',  ...
                       'FontSize', 10, ...
                       'String', 'Colormap',  ...
                       'Tag', 'color');
    enabled = [enabled hColor];

    hText = uicontrol('Parent', hPanel, ...
                      'Units', 'normalized',  ...
                      'Position', [0.89 0.69 0.075 0.05], ...
                      'String', 'Compression',  ...
                      'Style', 'text',  ...
                      'Tag', 'text19');

    hCompress = uicontrol('Parent', hPanel, ...
                          'Units', 'normalized',  ...
                          'Callback', @gui_Callback, ...
                          'Position', [0.89 0.655 0.075 0.05], ...
                          'String', typecompress, ...
                          'Style', 'popupmenu',  ...
                          'Value', 1, ...
                          'Tag', 'compression');
    enabled = [enabled hCompress];

    % The various filters, along with their labels
    hText = uicontrol('Parent', hPanel, ...
                      'Units', 'normalized',  ...
                      'Position', [0.9 0.525 0.05 0.05], ...
                      'String', 'Filters:',  ...
                      'FontSize', 12, ...
                      'FontWeight', 'bold', ...
                      'Style', 'text',  ...
                      'Tag', 'text16');

    hDetrend = uicontrol('Parent', hPanel, ...
                         'Units', 'normalized',  ...
                         'Callback', @gui_Callback, ...
                         'Position', [0.9 0.5 0.1 0.05], ...
                         'String', 'Detrend',  ...
                         'Style', 'checkbox',  ...
                         'Tag', 'detrend');
    enabled = [enabled hDetrend];

    hCosmics = uicontrol('Parent', hPanel, ...
                           'Units', 'normalized',  ...
                           'Callback', @gui_Callback, ...
                           'Position', [0.9 0.45 0.1 0.05], ...
                           'String', 'Cosmic rays',  ...
                           'Style', 'checkbox',  ...
                           'Tag', 'cosmics');
    enabled = [enabled hCosmics];

    hHotPixels = uicontrol('Parent', hPanel, ...
                         'Units', 'normalized',  ...
                         'Callback', @gui_Callback, ...
                         'Position', [0.9 0.4 0.1 0.05], ...
                         'String', 'Hot pixels',  ...
                         'Style', 'checkbox',  ...
                         'Tag', 'hot_pixels');
    enabled = [enabled hHotPixels];

    hNorm = uicontrol('Parent', hPanel, ...
                           'Units', 'normalized',  ...
                           'Callback', @gui_Callback, ...
                           'Position', [0.9 0.35 0.1 0.05], ...
                           'String', 'Normalize',  ...
                           'Style', 'checkbox',  ...
                           'Tag', 'normalize');
    enabled = [enabled hNorm];

    % The Snapshot button
    hSnapshot = uicontrol('Parent', hFig, ...
                    'Units', 'normalized',  ...
                    'Callback', @options_Callback, ...
                    'Position', [0.01 0.93 0.05 0.05], ...
                    'String', 'Snapshot',  ...
                    'Tag', 'snapshot');
    enabled = [enabled hSnapshot];

    % The buttons which allows to edit, load and save parameters
    hEdit = uicontrol('Parent', hPanel, ...
                       'Units', 'normalized',  ...
                       'Callback', @options_Callback, ...
                       'Position', [0.89 0.3 0.08 0.04], ...
                       'Style', 'pushbutton',  ...
                       'FontSize', 10, ...
                       'String', 'Edit parameters',  ...
                       'Tag', 'edit');
    enabled = [enabled hEdit];

    hLoad = uicontrol('Parent', hPanel, ...
                       'Units', 'normalized',  ...
                       'Callback', @options_Callback, ...
                       'Position', [0.89 0.2 0.08 0.04], ...
                       'Style', 'pushbutton',  ...
                       'FontSize', 10, ...
                       'String', 'Load parameters',  ...
                       'Tag', 'load');
    enabled = [enabled hLoad];

    hSave = uicontrol('Parent', hPanel, ...
                       'Units', 'normalized',  ...
                       'Callback', @options_Callback, ...
                       'Position', [0.89 0.15 0.08 0.04], ...
                       'Style', 'pushbutton',  ...
                       'FontSize', 10, ...
                       'String', 'Save parameters',  ...
                       'Tag', 'save');
    enabled = [enabled hSave];

    % We store all the useful handles into a structure to easily retrieve them,
    % along with some indexes
    handles = struct('uipanel', hPanel, ...
                     'slider', hFrame, ...
                     'text', hIndex, ...
                     'detrend', hDetrend, ...
                     'list', hChannel, ...
                     'hot_pixels', hHotPixels, ...
                     'cosmics', hCosmics, ...
                     'normalize', hNorm, ...
                     'channel_color', hColor, ...
                     'channel_type', hType, ...
                     'compress', hCompress, ...
                     'axes', [hAxes hAxesNext], ...
                     'experiment', hName, ...
                     'all_buttons', enabled, ...
                     'img', -1, ...
                     'prev_frame', -1, ...
                     'frame', 1, ...
                     'display', [1 1], ...
                     'prev_channel', -1, ...
                     'current', 1);

    % Link both axes to keep the same information on both sides
    linkaxes(handles.axes);

    return;
  end
end
