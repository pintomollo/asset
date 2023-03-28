function [myrecording, opts, is_updated] = inspect_segmentation(myrecording, opts)
% INSPECT_SEGMENTATION displays a pop-up window for the user to manually inspect the
% segmentation that will be performed on the provided movie.
%
%   [MYRECORDING, OPTS] = INSPECT_SEGMENTATION(MYRECORDING,OPTS) displays the window
%   using the data contained in MYRECORDING and the parameter values from OPTS. It
%   them accordingly to the user's choice. MYRECORDING should be a structure as
%   defined by get_struct('myrecording')
%
%   [...] = INSPECT_SEGMENTATION() prompts the user to select a MYRECORDING containing
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
  segmentations = get_struct('segmentation', [1, nchannels]);

  % Dragzoom help message
  imghelp = regexp(help('dragzoom'), ...
             '([ ]+Normal mode:.*\S)\s+Mouse actions in 3D','tokens');
  imghelp = ['DRAGZOOM interactions (help dragzoom):\n\n', imghelp{1}{1}];

  % Create the GUI using segmentations
  [hFig, handles] = create_figure();

  % Allocate the various variables. This allows them to be "persistent" between
  % different calls to the callback functions.
  img = [];
  orig_img = [];
  spots = [];
  filt_spots = [];
  reconstr = [];
  reconstr_filt = [];
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
    % Store the segmentations
    myrecording.segmentations = segmentations;
    % And get the experiment name
    myrecording.experiment = get(handles.experiment, 'String');
    % And reset the other fields
    for i=1:nchannels
      myrecording.segmentations(i).detections = get_struct('detection',0);
    end
    myrecording.trackings = get_struct('tracking', 0);
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

      % The filters
      set(handles.detrend,'Value', segmentations(indx).detrend);
      set(handles.denoise,'Value', segmentations(indx).denoise);
      set(handles.filter_spots,'Value', segmentations(indx).filter_spots);

      % The type of segmentation
      set(handles.segmentation_type, 'Value',  segmentations(indx).type);

      % And setup the indexes correctly
      handles.prev_channel = indx;
      handles.prev_frame = -1;
    end

    % Get the type of segmentation currently used
    contents = get(handles.segmentation_type,'String');
    segment_type = contents{segmentations(indx).type};

    % Here we actually segment the image (and update some important displays)
    if (recompute)

      % Because it takes long, display it and block the GUI
      set(hFig, 'Name', 'Images Segmentation (Processing...)');
      set(handles.all_buttons, 'Enable', 'off');
      drawnow;
      refresh(hFig);

      % Update the frame index
      set(handles.text, 'String', ['Frame #' num2str(nimg)]);

      % Delete any previous value for the noise level
      noise = [];

      % Load the new image
      orig_img = double(load_data(channels(indx).fname, nimg));

      % Copy it to the working variable
      img = orig_img;

      % Update the index
      handles.prev_frame = nimg;

      % Detrend the image ?
      if (segmentations(indx).detrend)
        img = imdetrend(img, opts.segmenting.detrend_meshpoints);
      end

      % Denoise the image ?
      if (segmentations(indx).denoise)
        [img, noise] = imdenoise(img, opts.segmenting.denoise_remove_bkg, ...
                        opts.segmenting.denoise_func, opts.segmenting.denoise_size);
      end

      % Segment the image and estimate the corresponding locations
      spots = perform_step('segmentation', segment_type, img, opts, noise);
      spots = perform_step('estimation', segment_type, img, spots, opts);

      % Estimate the noise if not previously done
      if (isempty(noise))
        noise = estimate_noise(img);
      end

      % And filter the spots
      filt_spots = perform_step('filtering', segment_type, spots, opts, noise);

      % Finally, reconstrcut the image using the previous detection
      reconstr = perform_step('reconstructing', segment_type, orig_img, spots);
      reconstr_filt = perform_step('reconstructing', segment_type, orig_img, filt_spots);
    end

    % Decide which type of spots to display
    if (segmentations(indx).filter_spots)
      curr_spots = filt_spots;
      curr_rec = reconstr_filt;
    else
      curr_spots = real(spots);
      curr_rec = reconstr;
    end

    % Determine which image to display in the left panel
    switch handles.display(1)

      % The original image
      case 2
        img1 = orig_img;

      % The difference between the filtered and the original images
      case 3
        img1 = orig_img - img;

      % The filtered image
      otherwise
        img1 = img;
    end

    % Determine which image to display in the right panel
    switch handles.display(2)

      % The reconstructed image
      case 2
        img2 = curr_rec;

      % The difference between the filtered and the reconstructed images
      case 3
        img2 = img - curr_rec;

      % The filtered image
      otherwise
        img2 = img;
    end

    % Get the colors for the display
    spots_colors = colors.spots{color_index};

    % If we have already created the axes and the images, we can simply change their
    % content (i.e. CData)
    if (numel(handles.img) > 1 & all(ishandle(handles.img)))
      set(handles.img(1),'CData', img1);
      set(handles.img(2),'CData', img2);

      % And update the spots
      perform_step('plotting', segment_type, handles.data, curr_spots, spots_colors);
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

      % Now add the detected spots
      handles.data = perform_step('plotting', segment_type, handles.axes(2), curr_spots, spots_colors);

      % Drag and Zoom library from Evgeny Pr aka iroln
      dragzoom(handles.axes, 'on');
    end

    % And set the colormap
    colormap(hFig, colors.colormaps{color_index}());

    if (recompute)
      % Release the image
      set(hFig, 'Name', 'Images Segmentation');
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
        [opts.segmenting, recompute] = edit_options(opts.segmenting);

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

    % By default we recompute the filter
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
      case {'detrend','denoise'}
        segmentations(indx).(type) = logical(get(hObject, 'Value'));

      % But filtering does not require recomputing
      case 'filter_spots'
        segmentations(indx).(type) = logical(get(hObject, 'Value'));
        recompute = false;

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
      case 'type'
        segmentations(indx).(type) = get(hObject, 'Value');

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
  % the segmentations structure for its standard form before releasing the
  % GUI to exit it.

    % Create a copy of segmentations in case we need to cancel
    tmp_segmentations = segmentations;

    % We need ot loop over the segmentations
    nsegmentations = length(segmentations);

    % Get the available types of segmentations and compressions
    contents = get(handles.segmentation_type,'String');
    ntypes = length(contents);

    % We want to check if some segmentations have identical features
    denoise = logical(zeros(nsegmentations,1));
    detrend = logical(zeros(nsegmentations,1));
    types = logical(zeros(nsegmentations,ntypes));

    % Convert the indexes into strings and build a summary of the filters
    for i=1:nsegmentations
      denoise(i) = segmentations(i).denoise;
      detrend(i) = (segmentations(i).detrend && channels(i).detrend);
      types(i,segmentations(i).type) = true;
      tmp_segmentations(i).type = contents{segmentations(i).type};
    end

    % Some checks to make sure the user is aware of some potential issues
    ok = true;
    if (ok & any(detrend, 1))

      % This is slow process which have apparently already been applied
      answer = questdlg('Some channels will be detrended a second time, continue ?');
      ok = strcmp(answer,'Yes');
    end
    if (ok & any(sum(types(:,2:end), 1)>1))

      % Just in case, for later display
      answer = questdlg('Multiple channels will be segmented, continue ?');
      ok = strcmp(answer,'Yes');
    end

    % If everything is OK, release the GUI and quit
    if (ok)
      segmentations = tmp_segmentations;
      uiresume(hFig);
    end

    return
  end

  function [hFig, handles] = create_figure
  % This function actually creates the GUI, placing all the elements
  % and linking the calbacks.

    % The number of channels provided
    nchannels = length(channels);

    % Initialize the possible segmentations and their corresponding channels
    typestring = {'None','multiscale_gaussian_spots', 'rectangular_local_maxima'};
    typechannel = {'luminescence','fluorescence'};

    % Initialize the structure used for the interface
    liststring = '';
    for i = 1:nchannels

      % Build the displayed list
      liststring = [liststring channels(i).type num2str(i)];
      if (i < nchannels)
        liststring = [liststring '|'];
      end

      % Set the segmentation type
      type_test = ismember(typechannel, channels(i).type);
      if any(type_test)
        segmentations(i).type = find(type_test)+1;
      else
        segmentations(i).type = 1;
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
                  'Name', 'Images Segmentation',  ...
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
                    'Max', nframes-1, ...
                    'Min', 1, ...
                    'Style', 'slider', ...
                    'Tag', 'slider');
    enabled = [enabled hFrame];

    %%%%%%% Now the main panel

    % The panel itsel
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
                         'String', 'Filtered image', ...
                         'Tag', 'radio21');
    enabled = [enabled hControl];

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.4 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Reconstructed image', ...
                         'Tag', 'radio22');
    enabled = [enabled hControl];

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.7 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Difference', ...
                         'Tag', 'radio23');
    enabled = [enabled hControl];

    % The type of segmentation to perform, along with its label
    hText = uicontrol('Parent', hPanel, ...
                      'Units', 'normalized',  ...
                      'Position', [0.875 0.925 0.075 0.05], ...
                      'String', 'Segmentation:',  ...
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

    hDenoise = uicontrol('Parent', hPanel, ...
                         'Units', 'normalized',  ...
                         'Callback', @gui_Callback, ...
                         'Position', [0.9 0.45 0.1 0.05], ...
                         'String', 'Denoise',  ...
                         'Style', 'checkbox',  ...
                         'Tag', 'denoise');
    enabled = [enabled hDenoise];

    hFilter = uicontrol('Parent', hPanel, ...
                         'Units', 'normalized',  ...
                         'Callback', @gui_Callback, ...
                         'Position', [0.9 0.4 0.1 0.05], ...
                         'String', 'Filter spots',  ...
                         'Style', 'checkbox',  ...
                         'Tag', 'filter_spots');
    enabled = [enabled hFilter];

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
                     'denoise', hDenoise, ...
                     'filter_spots', hFilter, ...
                     'list', hChannel, ...
                     'segmentation_type', hType, ...
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
