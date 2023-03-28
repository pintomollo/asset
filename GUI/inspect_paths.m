function [myrecording, opts, is_updated] = inspect_paths(myrecording, opts)
% INSPECT_PATHS displays a pop-up window for the user to manually inspect the
% filtering of the tracks that will be performed on the provided movie.
%
%   [MYRECORDING, OPTS] = INSPECT_PATHS(MYRECORDING,OPTS) displays the window using
%   the data contained in MYRECORDING and the parameter values from OPTS. It updates
%   them accordingly to the user's choice. MYRECORDING should be a structure as
%   defined by get_struct('myrecording')
%
%   [...] = INSPECT_PATHS() prompts the user to select a MYRECORDING containing
%   Matlab file before opening the GUI.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 17.06.2014

  % Argument checking, need to know if we ask for a myrecording file or not.
  if (nargin ~= 2 | isempty(myrecording) | isempty(opts))

    % Fancy output
    disp('[Select a MAT file]');

    % Prompting the user for the MAT file
    [fname, dirpath] = uigetfile({'*.mat'}, ['Load a MAT file']);
    fname = fullfile(dirpath, fname);

    % Loading was cancelled
    if (isequal(dirpath, 0))
      myrecording = [];
      opts = [];
      is_updated = false;

      return;
    end

    % Load the matrix and check its content
    data = load(fname);

    % Not what we expected
    if (~isfield(data, 'myrecording') || ~isfield(data, 'opts'))
      disp(['Error: ' fname ' does not contain a valid tracking structure']);
      myrecording = [];
      opts = [];
      is_updated = false;

      return;

    % Extract the loaded data
    else
      myrecording = data.myrecording;
      opts = data.opts;
    end
  end

  % Store the original options
  orig_opts = opts;

  % Prepare some global variables
  channels = myrecording.channels;
  nchannels = length(channels);
  segmentations = myrecording.segmentations;
  trackings = myrecording.trackings;

  % Copy the tracking data at the proper place
  for i=1:nchannels
    nframes = size_data(channels(i));
    if (length(trackings(i).filtered) < nframes)
      trackings(i).filtered = get_struct('detection', [1 nframes]);
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
  orig_img = [];
  all_paths = [];
  paths = [];
  all_colors = [];
  filtered_colors = [];
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

  if (is_updated)
    % Store the channels
    myrecording.channels = channels;
    % Store the trackings
    myrecording.trackings = trackings;
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
      color_index = channels(indx).color(1);

      % Set the name of the current panel
      set(handles.uipanel,'Title', [channels(indx).type ' ' num2str(indx)]);

      % And put the GUI in blocked mode
      set(hFig, 'Name', 'Tracks Filtering (Processing...)');
      set(handles.all_buttons, 'Enable', 'off');
      recompute = true;

      % The filters
      set(handles.refine,'Value', trackings(indx).reestimate_spots);

      % The paths
      all_paths = reconstruct_tracks(trackings(indx).detections, true);
      all_colors = colorize_graph(all_paths, colors.paths{color_index}(length(all_paths)));

      % And setup the indexes correctly
      handles.prev_channel = indx;
      handles.prev_frame = -1;
    end

    % Get the type of segmentation currently used
    segment_type = segmentations(indx).type;

    % The slider
    set(handles.text, 'String', ['Frame #' num2str(nimg)]);

    % Try to avoid reloading frames as much as possible
    if (handles.prev_frame ~= nimg)
      orig_img = double(load_data(channels(indx).fname, nimg));
    end

    % Here we filter the tracks
    if (recompute)

      % Because it takes long, display it and block the GUI
      set(hFig, 'Name', 'Tracks Filtering (Processing...)');
      set(handles.all_buttons, 'Enable', 'off');
      drawnow;
      refresh(hFig);

      % Perform the actual filtering of the paths
      links = filter_tracking(trackings(indx).detections, opts.tracks_filtering.min_path_length, opts.tracks_filtering.max_zip_length,logical(opts.tracks_filtering.interpolate));
      paths = reconstruct_tracks(links, true);
      filtered_colors = colorize_graph(paths, colors.paths{color_index}(length(paths)));
    end

    % Extract from the paths the current spots
    spots = cellfun(@(x)(x(x(:,end-1)==nimg,:)), all_paths, 'UniformOutput', false);
    spots = cat(1,spots{:});

    % And the filtered ones as well
    spots_filt = cellfun(@(x)(x(x(:,end-1)==nimg,:)), paths, 'UniformOutput', false);
    spots_filt = cat(1,spots_filt{:});

    % Do we need to reestimate the spots ?
    if (trackings(indx).reestimate_spots)

      % Fancy display
      if (~recompute)
        set(hFig, 'Name', 'Tracks Filtering (Processing...)');
        set(handles.all_buttons, 'Enable', 'off');
        drawnow;
        refresh(hFig);
      end

      % Perform the actual estimation
      spots_filt(:,2:end-2) = reestimate_spots(spots_filt(:,2:end-2), orig_img, segmentations(indx), opts);

      if (~recompute)
        set(hFig, 'Name', 'Tracks Filtering');
        set(handles.all_buttons, 'Enable', 'on');
      end
    end

    % Update the index
    handles.prev_frame = nimg;

    % Determine which type of image to display in the left panel
    switch handles.display(1)

      % Only the original links pointing towards the current frame
      case 2
        if (~isempty(spots))
          divs = spots(:,1);
          spots1 = {spots(divs<0,2:end), spots(divs==0,2:end), spots(divs>0,2:end)};
        else
          spots1 = {[]};
        end
        links1 = cellfun(@(x)(x(abs(x(:,end-1)-nimg) < 2,:)), all_paths, 'UniformOutput', false);
        links1 = links1(~cellfun('isempty', links1));
        colors1 = colorize_graph(links1, colors.paths{color_index}(length(links1)));

      % All the original full tracks
      case 3
        if (~isempty(spots))
          divs = spots(:,1);
          spots1 = {spots(divs<0,2:end), spots(divs==0,2:end), spots(divs>0,2:end)};
        else
          spots1 = {[]};
        end
        links1 = all_paths;
        colors1 = all_colors;

      % No paths at all
      otherwise
        spots1 = {[]};
        links1 = {[]};
        colors1 = 'k';
    end

    % Determine which type of image to display in the right panel
    switch handles.display(2)

      % Only the filtered links pointing towards the current frame
      case 2
        if (~isempty(spots_filt))
          divs = spots_filt(:,1);
          spots2 = {spots_filt(divs<0,2:end), spots_filt(divs==0,2:end), spots_filt(divs>0,2:end)};
        else
          spots2 = {[]};
        end

        % Keep only the links pointing at the current frame
        links2 = cellfun(@(x)(x(abs(x(:,end-1)-nimg) < 2,:)), paths, 'UniformOutput', false);
        links2 = links2(~cellfun('isempty', links2));
        colors2 = colorize_graph(links2, colors.paths{color_index}(length(links2)));

      % All the filtered full tracks
      case 3
        if (~isempty(spots_filt))
          divs = spots_filt(:,1);
          spots2 = {spots_filt(divs<0,2:end), spots_filt(divs==0,2:end), spots_filt(divs>0,2:end)};
        else
          spots2 = {[]};
        end
        links2 = paths;
        colors2 = filtered_colors;

      % No paths at all
      otherwise
        spots2 = {[]};
        links2 = {[]};
        colors2 = 'k';
    end

    % Get the colors for the display
    divisions_colors = colors.status{color_index};

    % If we have already created the axes and the images, we can simply change their
    % content (i.e. CData)
    if (numel(handles.img) > 1 & all(ishandle(handles.img)))
      % The images
      set(handles.img(1),'CData', orig_img);
      set(handles.img(2),'CData', orig_img);

      % The paths
      plot_paths(handles.data(3), links1, colors1);
      plot_paths(handles.data(4), links2, colors2);

      % And the spots on top
      perform_step('plotting', segment_type, handles.data(1), spots1, divisions_colors, true);
      perform_step('plotting', segment_type, handles.data(2), spots2, divisions_colors, true);
    else

      % Otherwise, we create the two images in their respective axes
      handles.img = image(orig_img,'Parent', handles.axes(1),...
                        'CDataMapping', 'scaled',...
                        'Tag', 'image');
      handles.img(2) = image(orig_img,'Parent', handles.axes(2), ...
                        'CDataMapping', 'scaled',...
                        'Tag', 'image');

      % Hide the axes and prevent a distortion of the image due to stretching
      set(handles.axes,'Visible', 'off',  ...
                 'DataAspectRatio',  [1 1 1]);

      % Now add the links
      handles.data(3) = plot_paths(handles.axes(1), links1, colors1);
      handles.data(4) = plot_paths(handles.axes(2), links2, colors2);

      % And their detected spots
      handles.data(1) = perform_step('plotting', segment_type, handles.axes(1), spots1, divisions_colors);
      handles.data(2) = perform_step('plotting', segment_type, handles.axes(2), spots2, divisions_colors);

      % Drag and Zoom library from Evgeny Pr aka iroln
      dragzoom(handles.axes, 'on')
    end

    % And set the colormap
    colormap(hFig, colors.colormaps{color_index}());

    if (recompute)
      % Release the image
      set(hFig, 'Name', 'Tracks Filtering');
      set(handles.all_buttons, 'Enable', 'on');
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
        [opts.tracks_filtering, recompute] = edit_options(opts.tracks_filtering);

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
      case 'reestimate_spots'
        trackings(indx).(type) = logical(get(hObject, 'Value'));

      % The slider varies the frame index
      case 'slider'
        handles.frame = round(get(hObject, 'Value'));
        recompute = false;

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
      case 'type'
        trackings(indx).(type) = get(hObject, 'Value');

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
  % This function releases the GUI to exit it

    uiresume(hFig);

    return
  end

  function [hFig, handles] = create_figure
  % This function actually creates the GUI, placing all the elements
  % and linking the callbacks.

    % The number of channels provided
    nchannels = length(channels);

    % Initialize the structure used for the interface
    liststring = '';
    for i = 1:nchannels

      % Build the displayed list
      liststring = [liststring channels(i).type num2str(i)];
      if (i < nchannels)
        liststring = [liststring '|'];
      end
    end

    % Get the number of frames
    nframes = size_data(channels(1).fname);

    % And the experiment name
    exp_name = myrecording.experiment;

    % Create my own grayscale map for the image display
    mygray = [0:255]' / 255;
    mygray = [mygray mygray mygray];

    % We build a list of all buttons to easily block and release them
    enabled = [];

    % The main figure, cannot be rescaled, closed nor deleted
    hFig = figure('PaperUnits', 'centimeters',  ...
                  'CloseRequestFcn', @gui_Callback, ...
                  'Color',  [0.7 0.7 0.7], ...
                  'Colormap', mygray, ...
                  'MenuBar', 'none',  ...
                  'Name', 'Tracks Filtering',  ...
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
                    'Max', max(nframes,1.1), ...
                    'Min', 1, ...
                    'Style', 'slider', ...
                    'Tag', 'slider');
    enabled = [enabled hFrame];

    %%%%%%% Now the main panel

    % The panel itself
    hPanel = uipanel('Parent', hFig, ...
                     'Title', [channels(1).type '1'],  ...
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
                         'String', 'No path', ...
                         'Tag', 'radio11');
    enabled = [enabled hControl];

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.4 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Current links', ...
                         'Tag', 'radio12');
    enabled = [enabled hControl];

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.7 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Full paths', ...
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
                         'String', 'No path', ...
                         'Tag', 'radio21');
    enabled = [enabled hControl];

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.4 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Filtered links', ...
                         'Tag', 'radio22');
    enabled = [enabled hControl];

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.7 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Filtered paths', ...
                         'Tag', 'radio23');
    enabled = [enabled hControl];

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

    % The various options for the filtering of the paths
    hRefine = uicontrol('Parent', hPanel, ...
                         'Units', 'normalized',  ...
                         'Callback', @gui_Callback, ...
                         'Position', [0.9 0.5 0.1 0.05], ...
                         'String', 'Reestimate spots',  ...
                         'Style', 'checkbox',  ...
                         'Tag', 'reestimate_spots');
    enabled = [enabled hRefine];

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
                     'refine', hRefine, ... %'force', hForce, ...
                     'list', hChannel, ...
                     'axes', [hAxes hAxesNext], ...
                     'experiment', hName, ...
                     'all_buttons', enabled, ...
                     'img', -1, ...
                     'data', -1, ...
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
