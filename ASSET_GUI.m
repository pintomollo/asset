function [] = ASSET_GUI(myrecording, opts)
% ASSET_GUI displays the main window of this interactive segmentation and
% tracking platform in time-lapse recordings.
%
%   ASSET_GUI() displays an empty GUI for the user to load a recording interactively.
%
%   [MYRECORDING, OPTS] = ASSET_GUI() returns the results of the analysis in MYRECORDING
%   and the corresponding options in OPTS. MYRECORDING and OPTS will be structured as defined
%   get_struct('myrecording') and get_struct('options') respectively.
%
%   [MYRECORDING, OPTS] = ASSET_GUI(MYRECORDING, OPTS) displays the window using
%   the data contained in MYRECORDING and the parameter values from OPTS. It updates
%   them accordingly to the user's choice.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 02.07.2014

  % Argument checking, need to know if we create new structures or not
  if (nargin ~= 2 || isempty(myrecording) || isempty(opts))
    myrecording = get_struct('myrecording');
    opts = get_struct('options');
  else
    % We utilize this function to improve compatibility between versions of this
    % platform, fusing option structures if need be.
    opts = update_structure(opts, 'options');
  end

  % Don't stay inside the directory itself
  if(exist('./install_ASSET.m', 'file'))
    cd('..')
  end

  % Create the GUI
  [hFig, handles] = create_figure();

  % Store all the useful data in a single structure to be stored in the figure
  %data = struct('channels', myrecording.channels, ...
  data = struct('recording', myrecording, ...
                'options', opts, ...
                'hFig', hFig, ...
                'img', [], ...
                'orig_img', [], ...
                'img_next', [], ...
                'handles', handles);

  % Update its content
  data = update_display(data);

  % Store the data structure in the main figure to be retrieved by the callbacks
  guidata(hFig, data);

  % Display the figure
  set(hFig,'Visible', 'on');

  return;
end

function [hFig, handles] = create_figure()
% This function actually creates the GUI, placing all the elements
% and linking the callbacks.

  % The number of channels provided
  %nchannels = length(channels);

  % Initialize the possible types and compressions
  %typestring = {'luminescence';'brightfield'; 'dic'; 'fluorescence'};
  %typecompress = {'none', 'lzw', 'deflate', 'jpeg'};

  % Initialize the structure used for the interface
  %liststring = '';
  %for i = 1:nchannels

  %  % Build the displayed list
  %  liststring = [liststring channels(i).type num2str(i)];
  %  if (i < nchannels)
  %    liststring = [liststring '|'];
  %  end
  %end

  % Get the number of frames
  %if nchannels > 0
  %  nframes = size_data(channels(1).fname);
  %else
  %  nframes = NaN;
  %end

  nframes = 1;

  % Set the parameters for the sliders
  %if isfinite(nframes)
  %  slider_step = [1 10]/nframes;
  %  slider_max = max(nframes, 1.1);
  %  slider_min = 1;
  %else
    slider_step = [1 1];
    slider_max = 1.1;
    slider_min = 1;
  %end

  % Create a name for the experiment based on the filename
  %exp_name = channels(1).fname;
  %[junk, exp_name, junk] = fileparts(exp_name);
  %[junk, exp_name, junk] = fileparts(exp_name);
  %exp_name = regexprep(exp_name, ' ', '');

  % Create my own grayscale map for the image display
  mygray = [0:255]' / 255;
  mygray = [mygray mygray mygray];

  % We build a list of all buttons to easily block and release them
  enabled = [];

  % Retrieve the resolution of the screen and infer the position of the GUI
  resol = get(0, 'ScreenSize');
  img_pos = ceil([0.05*resol(3:4) 0.9*resol(3:4)]);

  % The main figure, cannot be rescaled, closed nor deleted
  hFig = figure('PaperUnits', 'centimeters',  ...
                'CloseRequestFcn', @gui_CloseRequestFcn, ...
                'Colormap', mygray, ...
                'Color',  [0.7 0.7 0.7], ...
                'MenuBar', 'none',  ...
                'Name', 'ASSET',  ...
                'NumberTitle', 'off',  ...
                'Units', 'pixels', ...
                'Position', img_pos, ...
                'DeleteFcn', @gui_Callback, ...
                'HandleVisibility', 'callback',  ...
                'Tag', 'channel_fig', ...
                'UserData', [], ...
                'Visible', 'off');

  %%%%%% Now the buttons around the main panel

  % The list of channels
  hChannel = uicontrol('Parent', hFig, ...
                       'Units', 'normalized',  ...
                       'Callback', @gui_Callback, ...
                       'Position', [0.01 0.11 0.1 0.79], ...
                       %'String', liststring, ...
                         'String', '', ...
                       'Style', 'listbox',  ...
                       'Value', 1, ...
                       'Tag', 'channels');
  enabled = [enabled hChannel];

  % The OK button
  %hOK = uicontrol('Parent', hFig, ...
  %                'Units', 'normalized',  ...
  %                'Callback', @channel_CloseRequestFcn, ...
  %                'Position', [0.70 0.02 0.18 0.05], ...
  %                'String', 'OK',  ...
  %                'Tag', 'pushbutton11');
  %enabled = [enabled hOK];

  % The Quit button
  hQuit = uicontrol('Parent', hFig, ...
                  'Units', 'normalized',  ...
                  'Callback', @gui_CloseRequestFcn, ...
                  'Position', [0.90 0.02 0.08 0.05], ...
                  'String', 'Quit',  ...
                  'Tag', 'pushbutton12');
  enabled = [enabled hQuit];

  % The Load/Save buttons
  hLoadExp = uicontrol('Parent', hFig, ...
                  'Units', 'normalized',  ...
                  'Callback', @experiment_Callback, ...
                  'Position', [0.81 0.93 0.06 0.05], ...
                  'String', 'Load',  ...
                  'Tag', 'load');
  enabled = [enabled hLoadExp];

  hSaveExp = uicontrol('Parent', hFig, ...
                  'Units', 'normalized',  ...
                  'Callback', @experiment_Callback, ...
                  'Position', [0.87 0.93 0.06 0.05], ...
                  'String', 'Save',  ...
                  'Tag', 'save');

  hExport = uicontrol('Parent', hFig, ...
                  'Units', 'normalized',  ...
                  'Callback', @experiment_Callback, ...
                  'Position', [0.93 0.93 0.06 0.05], ...
                  'String', 'Export',  ...
                  'Tag', 'export');

  % The experiment name and its labels
  hText = uicontrol('Parent', hFig, ...
                    'Units', 'normalized',  ...
                    'Position', [0.2 0.93 0.09 0.025], ...
                    'String', 'Experiment name:',  ...
                    'BackgroundColor', get(hFig, 'Color'), ...
                    'FontSize', 12, ...
                    'Style', 'text',  ...
                    'Tag', 'text1');

  hName = uicontrol('Parent', hFig, ...
                    'Units', 'normalized',  ...
                    'Position', [0.3 0.93 0.5 0.05], ...
                    %'String', exp_name,  ...
                      'String', 'Load a recording here -->', ...
                    'FontSize', 12, ...
                    'Style', 'edit',  ...
                    'Tag', 'experiment');
  enabled = [enabled hName];

  % The sliders and their labels
  hIndex1 = uicontrol('Parent', hFig, ...
                    'Units', 'normalized',  ...
                    'Position', [0.1 0.03 0.08 0.025], ...
                    'String', 'Frame #1',  ...
                    'BackgroundColor', get(hFig, 'Color'), ...
                    'FontSize', 12, ...
                    'Style', 'text',  ...
                    'Tag', 'text1');

  hFrame1 = uicontrol('Parent', hFig, ...
                    'Units', 'normalized',  ...
                    'Callback', @gui_Callback, ...
                    'Position', [0.19 0.03 0.28 0.025], ...
                    'Value', 1, ...
                    'SliderStep', [1 10]/nframes, ...
                    'Max', max(nframes, 1.1), ...
                    'Min', 1, ...
                    'Style', 'slider', ...
                    'Tag', 'slider1');
  enabled = [enabled hFrame1];

  hIndex2 = uicontrol('Parent', hFig, ...
                    'Units', 'normalized',  ...
                    'Position', [0.49 0.03 0.08 0.025], ...
                    'String', 'Frame #1',  ...
                    'BackgroundColor', get(hFig, 'Color'), ...
                    'FontSize', 12, ...
                    'Style', 'text',  ...
                    'Tag', 'text2');

  hFrame2 = uicontrol('Parent', hFig, ...
                  'Units', 'normalized',  ...
                  'Callback', @gui_Callback, ...
                  'Position', [0.58 0.03 0.28 0.025], ...
                  'Value', 1, ...
                  'SliderStep', [1 10]/nframes, ...
                  'Max', max(nframes, 1.1), ...
                  'Min', 1, ...
                  'Style', 'slider', ...
                  'Tag', 'slider2');
  enabled = [enabled hFrame2];

  %% The Add and Remove buttons
  %hAdd = uicontrol('Parent', hFig, ...
  %                'Units', 'normalized',  ...
  %                'Callback', @add_channel_Callback, ...
  %                'Position', [0.01 0.055 0.1 0.04], ...
  %                'String', 'Add channel',  ...
  %                'Tag', 'pushbutton13');
  %enabled = [enabled hAdd];

  %hRemove = uicontrol('Parent', hFig, ...
  %                'Units', 'normalized',  ...
  %                'Callback', @remove_channel_Callback, ...
  %                'Position', [0.01 0.01 0.1 0.04], ...
  %                'String', 'Remove channel',  ...
  %                'Tag', 'pushbutton14');
  %enabled = [enabled hRemove];

  %%%%%%% Now the main panel

  % The panel itsel
  hPanel = uipanel('Parent', hFig, ...
                   'Title', 'Image 1',  ...
                   'Tag', 'uipanel',  ...
                   'Clipping', 'on',  ...
                   'Position', [0.12 0.11 0.87 0.8]);

  % The two axes
  hAxes = axes('Parent', hPanel, ...
               'Position', [0.01 0.1 0.43 0.9], ...
               'DataAspectRatio', [1 1 1], ...
               'Visible', 'off',  ...
               'Tag', 'axes');

  hAxesNext = axes('Parent', hPanel, ...
               'Position', [0.45 0.1 0.43 0.9], ...
               'DataAspectRatio', [1 1 1], ...
               'Visible', 'off',  ...
               'Tag', 'axes');

    % The two radio button groups that handle which image to display
    % For the choices to be mutually exclusive, one has to put them inside
    % such uibuttongroup.
    hRadio = uibuttongroup('Parent', hPanel, ...
                         'Units', 'normalized',  ...
                         'SelectionChangedFcn', @gui_Callback, ...
                         'Position', [0 0.05 0.43 0.05], ...
                         'tag', 'radio');

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.02 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Image', ...
                         'Tag', 'radio11');
    enabled = [enabled hControl];

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.27 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Detections', ...
                         'Tag', 'radio12');
    enabled = [enabled hControl];

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.53 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Links', ...
                         'Tag', 'radio13');
    enabled = [enabled hControl];

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.78 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Paths', ...
                         'Tag', 'radio14');
    enabled = [enabled hControl];

    % The second group
    hRadio = uibuttongroup('Parent', hPanel, ...
                         'Units', 'normalized',  ...
                         'SelectionChangedFcn', @gui_Callback, ...
                         'Position', [0.44 0.05 0.43 0.05], ...
                         'tag', 'radio');

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.02 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Image', ...
                         'Tag', 'radio21');
    enabled = [enabled hControl];

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.27 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Detections', ...
                         'Tag', 'radio22');
    enabled = [enabled hControl];

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.53 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Links', ...
                         'Tag', 'radio23');
    enabled = [enabled hControl];

    hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'normalized',  ...
                         'Position', [0.78 0.1 0.25 0.8], ...
                         'Style', 'radiobutton',  ...
                         'String', 'Paths', ...
                         'Tag', 'radio24');
    enabled = [enabled hControl];

    % The New/Filter/Segment/Track buttons
    hNew = uicontrol('Parent', hPanel, ...
                    'Units', 'normalized',  ...
                    'Callback', @pipeline_Callback, ...
                    'Position', [0.89 0.88 0.09 0.06], ...
                    'String', 'New Experiment',  ...
                    'Tag', 'new');
    enabled = [enabled hNew];

    hProcess = uicontrol('Parent', hPanel, ...
                    'Units', 'normalized',  ...
                    'Callback', @pipeline_Callback, ...
                    'Position', [0.89 0.78 0.09 0.06], ...
                    'String', 'Load Recording',  ...
                    'Tag', 'process');
    enabled = [enabled hProcess];

    hSegment = uicontrol('Parent', hPanel, ...
                    'Units', 'normalized',  ...
                    'Callback', @pipeline_Callback, ...
                    'Position', [0.89 0.72 0.09 0.06], ...
                    'String', 'Segment eggshell',  ...
                    'Enable', 'off', ...
                    'Visible', 'off', ...
                    'Tag', 'segment');

    hTrack = uicontrol('Parent', hPanel, ...
                    'Units', 'normalized',  ...
                    'Callback', @pipeline_Callback, ...
                    'Position', [0.89 0.66 0.09 0.06], ...
                    'String', 'Segment membrane',  ...
                    'Enable', 'off', ...
                    'Visible', 'off', ...
                    'Tag', 'track');

    hPaths = uicontrol('Parent', hPanel, ...
                    'Units', 'normalized',  ...
                    'Callback', @pipeline_Callback, ...
                    'Position', [0.89 0.60 0.09 0.06], ...
                    'String', 'Quantification',  ...
                    'Enable', 'off', ...
                    'Visible', 'off', ...
                    'Tag', 'paths');

    hAutosave = uicontrol('Parent', hPanel, ...
                         'Units', 'normalized',  ...
                         'Callback', @gui_Callback, ...
                         'Position', [0.9 0.54 0.1 0.05], ...
                         'String', 'Autosave',  ...
                         'Style', 'checkbox',  ...
                         'Value', 1, ...
                         'Tag', 'autosave');
    enabled = [enabled hAutosave];

    % The buttons which allows to change the colormap
    hColor = uicontrol('Parent', hPanel, ...
                       'Units', 'normalized',  ...
                       'Callback', @gui_Callback, ...
                       'Position', [0.89 0.4 0.08 0.04], ...
                       'Style', 'pushbutton',  ...
                       'FontSize', 10, ...
                       'String', 'Colormap',  ...
                       'Tag', 'color');
    enabled = [enabled hColor];


  % The Add and Remove buttons
  %hControl = uicontrol('Parent', hPanel, ...
  %                'Units', 'normalized',  ...
  %                'Callback', @add_ROI_Callback, ...
  %                'Position', [0.125 0.01 0.2 0.075], ...
  %                'String', 'Add ROI',  ...
  %                'Tag', 'addroi');
  %enabled = [enabled hControl];

  % The Add and Remove buttons
  %hControl = uicontrol('Parent', hPanel, ...
  %                'Units', 'normalized',  ...
  %                %'Callback', @find_zooids_Callback, ...
  %                'Position', [0.565 0.01 0.2 0.075], ...
  %                'String', 'Find zooids',  ...
  %                'Tag', 'addroi');
  %enabled = [enabled hControl];

  %% The type, color and compression of the channel, along with its labels
  %hText = uicontrol('Parent', hPanel, ...
  %                  'Units', 'normalized',  ...
  %                  'Position', [0.9 0.925 0.05 0.05], ...
  %                  'String', 'Channel:',  ...
  %                  'FontSize', 12, ...
  %                  'FontWeight', 'bold', ...
  %                  'Style', 'text',  ...
  %                  'Tag', 'text16');

  %hText = uicontrol('Parent', hPanel, ...
  %                  'Units', 'normalized',  ...
  %                  'Position', [0.885 0.875 0.075 0.05], ...
  %                  'String', 'Pixel size',  ...
  %                  'TooltipString', 'Pixel resolution in um', ...
  %                  'Style', 'text',  ...
  %                  'Tag', 'text17');

  %hResol = uicontrol('Parent', hPanel, ...
  %                  'Units', 'normalized',  ...
  %                  'Position', [0.89 0.845 0.07 0.04], ...
  %                  'String', '2.5', ...
  %                  'Style', 'edit',  ...
  %                  'Tag', 'data');
  %enabled = [enabled hResol];

  %hText = uicontrol('Parent', hPanel, ...
  %                  'Units', 'normalized',  ...
  %                  'Position', [0.9 0.77 0.05 0.05], ...
  %                  'String', 'Amplitude',  ...
  %                  'TooltipString', 'Intensity threshold required between two zooids (-1 = automatic threshold)', ...
  %                  'Style', 'text',  ...
  %                  'Tag', 'text17');

  %hAmpl = uicontrol('Parent', hPanel, ...
  %                  'Units', 'normalized',  ...
  %                  'Position', [0.89 0.74 0.07 0.04], ...
  %                  'Style', 'edit',  ...
  %                  'String', '-1', ...
  %                  'Tag', 'data');
  %enabled = [enabled hAmpl];

  %% The various filters, along with their labels
  %hText = uicontrol('Parent', hPanel, ...
  %                  'Units', 'normalized',  ...
  %                  'Position', [0.9 0.425 0.05 0.05], ...
  %                  'String', 'Filters:',  ...
  %                  'FontSize', 12, ...
  %                  'FontWeight', 'bold', ...
  %                  'Style', 'text',  ...
  %                  'Tag', 'text16');

  %hNorm = uicontrol('Parent', hPanel, ...
  %                   'Units', 'normalized',  ...
  %                   'Callback', @gui_Callback, ...
  %                   'Position', [0.9 0.35 0.1 0.05], ...
  %                   'String', 'Normalize',  ...
  %                   'Style', 'checkbox',  ...
  %                   'Tag', 'normalize');
  %enabled = [enabled hNorm];

  % The Snapshot button
  hSnapshot = uicontrol('Parent', hFig, ...
                  'Units', 'normalized',  ...
                  'Callback', @options_Callback, ...
                  'Position', [0.01 0.93 0.05 0.05], ...
                  'String', 'Snapshot',  ...
                  'Tag', 'snapshot');
  enabled = [enabled hSnapshot];

  % The Help button
  hHelp = uicontrol('Parent', hFig, ...
                  'Units', 'normalized',  ...
                  'Callback', @options_Callback, ...
                  'Position', [0.07 0.93 0.05 0.05], ...
                  'String', 'Help',  ...
                  'Tag', 'help');
  enabled = [enabled hHelp];

    % The Cleaning button
    hClean = uicontrol('Parent', hFig, ...
                    'Units', 'normalized',  ...
                    'Callback', @options_Callback, ...
                    'Position', [0.01 0.03 0.075 0.035], ...
                    'String', 'Clean TmpData',  ...
                    'Tag', 'clean');
    enabled = [enabled hClean];

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
                   'list', hChannel, ...
                   %'normalize', hNorm, ...
                   'axes', [hAxes hAxesNext], ...
                  'sliders', [hIndex1, hIndex2], ...
                   'autosave', true, ...
                   'experiment', hName, ...
                   'all_buttons', enabled, ...
                   %'resolution', hResol, ...
                   %'amplitude', hAmpl, ...
                   'hFig', hFig, ...
                    'img', -1, ...
                   'scatter', -1, ...
                   'roi', {{}}, ...
                   'display', [1 1], ...
                   'prev_channel', -1, ...
                    'prev_frame', [-1 -1], ...
                    'frame', [1 1], ...
                    'prev_channel', -1, ...
                   'current', -1);

  % Link both axes to keep the same information on both sides
  %linkaxes(handles.axes);

  return;
end

function [data] = handle_rois(data, prev_indx, next_indx)

  handles = data.handles;

  % Extract the ROIs
  %nrois = length(handles.roi);
  %if (prev_indx > 0)
  %  systems = NaN(0, 2);
  %  for i=1:nrois
  %    systems = [systems; handles.roi{i}.getPosition(); NaN(1,2)];
  %  end
  %  data.recording.channels(prev_indx).system = systems;
  %end

  %if (next_indx > 0)

  %  systems = data.recording.channels(next_indx).system;
  %  sys = find(all(isnan(systems), 2));
  %  nsys = length(sys);

  %  prev = 1;
  %  for i=1:nsys
  %    curr_sys = systems(prev:sys(i)-1,:);
  %    if (i > nrois)
  %      handles.roi{i} = impoly(handles.axes(1), curr_sys, 'Closed', false);
  %    else
  %      handles.roi{i}.setPosition(curr_sys);
  %    end
  %    prev = sys(i)+1;
  %  end

  %  for i=nrois:-1:nsys+1
  %    handles.roi{i}.delete();
  %    handles.roi(i) = [];
  %  end
  %end

  %data.handles = handles;

  return;
end

function data = update_display(data, recompute)
% The main figure of the GUI, the one responsible for the proper display
% of its content.

  % By default we recompute everything
  if (nargin < 2)
    recompute = true;
  end

  % Retrive the stored data
  handles = data.handles;
  channels = data.recording.channels;
  segmentations = data.recording.segmentations;
  trackings = data.recording.trackings;
  opts = data.options;
  hFig = handles.hFig;
  img = data.img;
  orig_img = data.orig_img;
  img_next = data.img_next;
  colors = get_struct('colors');

  % Get the indexes of the current frame and channel
  indx = handles.current;
  nimg = handles.frame;
  nchannels = length(channels);

  % Stop if no data at all
  if (indx < 1 || indx > nchannels)

    if (numel(handles.img) > 1 && all(ishandle(handles.img)))
      set(handles.img(1),'CData', []);
      set(handles.img(2),'CData', []);

      delete(get(handles.data(1), 'Children'));
      delete(get(handles.data(2), 'Children'));
      delete(get(handles.data(3), 'Children'));
      delete(get(handles.data(4), 'Children'));

      set(handles.scale, 'XData', [], 'YData', []);
    end

    return;
  end

  % If we have changed channel, we need to update the display of the buttons
  if (indx ~= handles.prev_channel && indx <= nchannels)
    % Because it takes long, display it and block the GUI
    color_index = channels(indx).color(1);
    set(hFig, 'Name', 'ASSET (Processing...)');
    %all_all = [handles.all_buttons, handles.save, handles.pipeline];
    %curr_status = get(all_all, 'Enable');
    set(handles.all_buttons, 'Enable', 'off');
    drawnow;
    refresh(hFig);

    has_segmentation = false;

    % The name
    typestring = {'luminescence';'brightfield'; 'dic'; 'fluorescence'};
    set(handles.uipanel,'Title', [typestring{channels(indx).type} ' ' num2str(indx)]);

    % And setup the indexes correctly
    handles.prev_channel = indx;
    handles.prev_frame = [-1 -1];

    % The paths
    all_paths = [];

    % Check what is available in the structure
    has_segmentation = false;
    has_tracking = false;
    has_filtered = false;

    if (~isempty(segmentations))
      has_segmentation = true;
      if (~isempty(trackings) && (length(trackings(indx).detections)==nframes))
        has_tracking = true;

        if (isfield(trackings(indx), 'filtered') && (length(trackings(indx).filtered)==nframes))
          for i=1:nframes
            has_filtered = (~isempty(trackings(indx).filtered(i).cluster));
            if has_filtered
              break;
            end
          end
        end
      end
    end

    % Reconstruct the available tracks
    if ~has_tracking
      all_paths = {[]};
    else
      all_paths = reconstruct_tracks(trackings(indx).detections, true);
    end
    if ~has_filtered
      all_filtered = {[]};
    else
      all_filtered = reconstruct_tracks(trackings(indx).filtered, true);
    end

    set(hFig, 'Name', 'ASSET');
    %set(handles.all_buttons, {'Enable'}, curr_status);
  end

  % And get the corresponding colormaps
  all_colors = colorize_graph(all_paths, colors.paths{color_index}(length(all_paths)));
  all_colors_filtered = colorize_graph(all_filtered, colors.paths{color_index}(length(all_filtered)));

  % The slider
  set(handles.sliders(1), 'String', ['Frame #' num2str(nimg(1))]);
  set(handles.sliders(2), 'String', ['Frame #' num2str(nimg(2))]);

  if (recompute)
    % Because it takes long, display it and block the GUI
    set(hFig, 'Name', 'ASSET (Processing...)');
    %curr_status = get(handles.all_buttons, 'Enable');
    set(handles.all_buttons, 'Enable', 'off');
    drawnow;
    refresh(hFig);

    % Here we recompute all the filtering of the frame
    noise = [];

    % Try to avoid reloading frames as much as possible
    if (nimg(1) == handles.prev_frame(1))
      if (nimg(2) == nimg(1))
        img_next = orig_img;
      else
        img_next = double(load_data(channels(indx).fname, nimg(2)));
      end
    elseif (nimg(2) == handles.prev_frame(2))
      if (nimg(2) == nimg(1))
        orig_img = img_next;
      else
        orig_img = double(load_data(channels(indx).fname, nimg(1)));
      end
    else
      if (nimg(2) == nimg(1))
        orig_img = double(load_data(channels(indx).fname, nimg(1)));
        img_next = orig_img;
      else
        orig_img = double(load_data(channels(indx).fname, nimg(1)));
        img_next = double(load_data(channels(indx).fname, nimg(2)));
      end
    end

    % Update the index
    handles.prev_frame = nimg;
  end

  % Determine which data to display in the left panel
  switch handles.display(1)

    % The segmented image
    case 2
      if has_segmentation
        spots1 = [segmentations(indx).detections(nimg(1)).carth segmentations(indx).detections(nimg(1)).properties];
      else
        spots1 = [];
      end
      links1 = {[]};
      colors1 = [];
      divisions_colors1 = colors.spots{color_index};

    % The tracked spots
    case 3

      if has_tracking
        spots = cellfun(@(x)(x(x(:,end-1)==nimg(1),:)), all_paths, 'UniformOutput', false);
        spots = cat(1,spots{:});
        divs = spots(:,1);
        spots1 = {spots(divs<0,2:end), spots(divs==0,2:end), spots(divs>0,2:end)};
        links1 = cellfun(@(x)(x(abs(x(:,end-1)-nimg(1)) < 2,:)), all_paths, 'UniformOutput', false);
        links1 = links1(~cellfun('isempty', links1));
      else
        spots1 = {[]};
        links1 = {[]};
      end
      colors1 = colorize_graph(links1, colors.paths{color_index}(length(links1)));
      divisions_colors1 = colors.status{color_index};

    % The full paths
    case 4
      if has_filtered
        spots = cellfun(@(x)(x(x(:,end-1)==nimg(1),:)), all_filtered, 'UniformOutput', false);
        spots = cat(1,spots{:});
        divs = spots(:,1);
        spots1 = {spots(divs<0,2:end), spots(divs==0,2:end), spots(divs>0,2:end)};
      else
        spots1 = {[]};
      end
      links1 = all_filtered;
      colors1 = all_colors_filtered;
      divisions_colors1 = colors.status{color_index};

    % The image only
    otherwise
      spots1 = {[]};
      links1 = {[]};
      colors1 = 'k';
      divisions_colors1 = 'k';
  end

  % Determine which image to display in the right panel
  switch handles.display(2)

    % The segmented image
    case 2
      if has_segmentation
        spots2 = [segmentations(indx).detections(nimg(2)).carth segmentations(indx).detections(nimg(2)).properties];
      else
        spots2 = {[]};
      end
      links2 = {[]};
      colors2 = [];
      divisions_colors2 = colors.spots_next{color_index};

    % The tracked data
    case 3
      if has_tracking
        spots_next = cellfun(@(x)(x(x(:,end-1)==nimg(2),:)), all_paths, 'UniformOutput', false);
        spots_next = cat(1,spots_next{:});
        divs = spots_next(:,1);
        spots2 = {spots_next(divs<0,2:end), spots_next(divs==0,2:end), spots_next(divs>0,2:end)};
        links2 = cellfun(@(x)(x(abs(x(:,end-1)-nimg(2)) < 2,:)), all_paths, 'UniformOutput', false);
        links2 = links2(~cellfun('isempty', links2));
      else
        spots2 = {[]};
        links2 = {[]};
      end
      colors2 = colorize_graph(links2, colors.paths{color_index}(length(links2)));
      divisions_colors2 = colors.status{color_index};

    % The full paths
    case 4
      if has_filtered
        spots_next = cellfun(@(x)(x(x(:,end-1)==nimg(2),:)), all_filtered, 'UniformOutput', false);
        spots_next = cat(1,spots_next{:});
        divs = spots_next(:,1);
        spots2 = {spots_next(divs<0,2:end), spots_next(divs==0,2:end), spots_next(divs>0,2:end)};
      else
        spots2 = {[]};
      end
      links2 = all_filtered;
      colors2 = all_colors;
      divisions_colors2 = colors.status{color_index};

    % The image alone
    otherwise
      spots2 = {[]};
      links2 = {[]};
      colors2 = 'k';
      divisions_colors2 = 'k';
  end

  % Get the type of segmentation used, if one is available
  if (has_segmentation)
    segment_type = segmentations(indx).type;
  else
    segment_type = 'unknown';
  end

  % If we have already created the axes and the images, we can simply change their
  % content (i.e. CData, XData, ...)
  [size_y, size_x] = size(orig_img);
  if (numel(handles.img) > 1 && all(ishandle(handles.img)))
    set(handles.img(1),'CData', orig_img);
    set(handles.img(2),'CData', img_next);

    plot_paths(handles.data(3), links1, colors1);
    plot_paths(handles.data(4), links2, colors2);

    perform_step('plotting', segment_type, handles.data(1), spots1, divisions_colors1, iscell(spots1));
    perform_step('plotting', segment_type, handles.data(2), spots2, divisions_colors2, iscell(spots2));

    set(handles.scale, 'XData', size_x*[0.05 0.05]+[0 10/opts.pixel_size], 'YData', size_y*[0.95 0.95]);
  else

    % Otherwise, we create the two images in their respective axes
    handles.img = image(orig_img,'Parent', handles.axes(1),...
                      'CDataMapping', 'scaled',...
                      'Tag', 'image');
    handles.img(2) = image(img_next,'Parent', handles.axes(2), ...
                      'CDataMapping', 'scaled',...
                      'Tag', 'image');

    % Hide the axes and prevent a distortion of the image due to stretching
    set(handles.axes,'Visible', 'off',  ...
               'DataAspectRatio',  [1 1 1]);

    % Now add the links
    handles.data(3) = plot_paths(handles.axes(1), links1, colors1);
    handles.data(4) = plot_paths(handles.axes(2), links2, colors2);

    % And their detected spots
    handles.data(1) = perform_step('plotting', segment_type, handles.axes(1), spots1, divisions_colors1, iscell(spots1));
    handles.data(2) = perform_step('plotting', segment_type, handles.axes(2), spots2, divisions_colors2, iscell(spots2));

    % And the necessary scale bar
    handles.scale = line('XData', size_x*[0.05 0.05]+[0 10/opts.pixel_size], 'YData', size_y*[0.95 0.95], 'Parent', handles.axes(1), 'Color', 'w', 'LineWidth', 4);

    % Drag and Zoom library from Evgeny Pr aka iroln
    dragzoom2D(handles.axes(1), 'on')
  end
  colormap(hFig, colors.colormaps{color_index});

  if (recompute || indx ~= handles.prev_channel)
    % Release the image
    set(hFig, 'Name', 'ASSET');
    set(handles.all_buttons, 'Enable', 'on');
  end

  % Update the data
  data.handles = handles;

  return
end

function remove_channel_Callback(hObject, eventdata)
% This function removes ones channel from the list

  data = guidata(hObject);
  handles = data.handles;

  % Get the current index for later
  indx = handles.current;

  % And ask for confirmation
  answer = questdlg(['Are you sure you want to remove channel ' num2str(indx) ' ?' ...
              '(No data will be deleted from the disk)'], 'Removing a channel ?');
  ok = strcmp(answer, 'Yes');

  % If it's ok, let's go
  if (ok)

    % Remove the current index and select the first one
    data.recording.channels(indx) = [];
    handles.current = 1;

    % If it was the only one, we need to handle this
    if (isempty(data.recording.channels))

      % Set up the indexes as empty
      handles.current = 0;
      handles.prev_channel = -1;

      % As this is long, block the GUI
      set(handles.all_buttons, 'Enable', 'off');
      set(handles.img, 'CData', []);
      set(handles.list, 'String', '');

      % And call the movie conversion function to load a new one
      new_channel = load_images();

      % Did we get something back ?
      if (~isempty(new_channel))

        % We provide basic default values for all fields
        data.recording.channels = get_struct('channel');
        data.recording.channels.fname = new_channel;

        % We update the list of available channels
        set(handles.list, 'String', 'Image 1', 'Value', 1);
      end

    % Otherwise, delete the last traces of the deleted channel
    else
      handles.prev_channel = -1;
      tmp_list = get(handles.list, 'String');
      tmp_list(indx,:) = [];
      set(handles.list, 'String', tmp_list, 'Value', 1);
    end

    % Update the data
    data.handles = handles;

    % Release the GUI and update
    data = update_display(data, true);

    guidata(handles.hFig, data);

    set(handles.all_buttons, 'Enable', 'on');
  end

  return;
end

function add_ROI_Callback(hObject, eventdata)

  data = guidata(hObject);
  handles = data.handles;

  % As this is long, block the GUI
  set(handles.all_buttons, 'Enable', 'off');
  set(handles.hFig, 'Name', 'ASSET (Draw your ROI...)');

  h = impoly(handles.axes(1), 'Closed', false);
  if (~isempty(h))
    handles.roi{end+1} = h;
  end

  % Update the data
  data.handles = handles;
  guidata(handles.hFig, data);

  % Release the GUI
  set(handles.hFig, 'Name', 'ASSET');
  set(handles.all_buttons, 'Enable', 'on');
end

%function find_zooids_Callback(hObject, eventdata)

%  data = guidata(hObject);
%  handles = data.handles;

  % As this is long, block the GUI
%  set(handles.all_buttons, 'Enable', 'off');

%  params = [8/str2double(get(handles.resolution, 'String')) str2double(get(handles.amplitude, 'String'))];

%  data = handle_rois(data, handles.current, -1);

  %keyboard
%
%  zooids = find_zooids(data.img, data.recording.channels(handles.current).system, params);
%  data.recording.channels(handles.current).zooids = zooids;


%  % Release the GUI
%  data = update_display(data, false);
%  guidata(handles.hFig, data);
%  set(handles.all_buttons, 'Enable', 'on');
%end

function add_channel_Callback(hObject, eventdata)
% This function adds a new channel to the current recording

  data = guidata(hObject);
  handles = data.handles;

  % Remember how many there were
  nchannels = length(data.recording.channels);

  % As this is long, block the GUI
  set(handles.all_buttons, 'Enable', 'off');

  % And call the movie conversion function to get a new channel
  new_channel = load_images();

  % Did we get anything ?
  if (~isempty(new_channel))

    liststring = get(handles.list, 'String');
    for i=1:length(new_channel)
      % We provide basic default values for all fields
      data.recording.channels(end+1) = get_struct('channel');
      data.recording.channels(end).fname = new_channel{i};

      % We update the list of available channels
      if (nchannels == 0)
        liststring = ['Image ' num2str(length(data.recording.channels))];
      else
        if (size(liststring, 1) > 1)
          liststring(:,end+1) = '|';
          liststring = liststring.';
          liststring = liststring(:).';
          liststring = liststring(1:end-1);
        end
        liststring = [liststring '|Image ' num2str(length(data.recording.channels))];
      end

      nchannels = nchannels + 1;
    end
    set(handles.list, 'String', liststring);
  end

  % If we just add a new first channel, we need to set it up properly and update
  if (nchannels == 0 && length(channels) > 0)
    handles.current = 1;
    set(handles.list, 'Value', 1);
    data.handles = handles;
    data = update_display(data, true);
  end

  % Release the GUI
  guidata(handles.hFig, data);
  set(handles.all_buttons, 'Enable', 'on');

  return
end

function pipeline_Callback(hObject, eventdata)
% This function is responsible for handling the content of the
% structure which contains the parameters of the filtering algorithms.

  data = guidata(hObject);
  handles = data.handles;
  myrecording = data.recording;
  opts = data.options;

  % Block the GUI
  set(handles.all_buttons, 'Enable', 'off');
  drawnow;
  refresh(handles.hFig);

  % And get the type of button which called the callback (from its tag)
  type = get(hObject, 'tag');

  % By default, recompute
  reload = true;

  % Handle all three buttons differently
  switch type

    % Create a new, empty experiment
    case 'new'
      myrecording = get_struct('myrecording');
      opts = get_struct('options');

    % Call the data processing GUI and process the channels accordingly
    case 'process'
      set(handles.hFig, 'Visible', 'off')
      [myrecording, opts, reload] = inspect_recording(myrecording, opts);
      if (reload)
        [myrecording, opts] = preprocess_movie(myrecording, opts);
        if (handles.autosave)
          save([myrecording.experiment '.mat'], 'myrecording', 'opts');
        end
        [opts, recompute] = edit_options(opts);
      end
      set(handles.hFig, 'Visible', 'on')

    % Call the segmenting GUI, and segment accordingly
    case 'segment'
      set(handles.hFig, 'Visible', 'off')
      [myrecording, opts, reload] = inspect_segmentation(myrecording, opts);
      if (reload)
        [myrecording, opts] = segment_movie(myrecording, opts);
        if (handles.autosave)
          save([myrecording.experiment '.mat'], 'myrecording', 'opts');
        end
      end
      set(handles.hFig, 'Visible', 'on')

    % Call the cell tracking GUI, and track accordingly
    case 'track'
      set(handles.hFig, 'Visible', 'off')
      [myrecording, opts, reload] = inspect_tracking(myrecording, opts);
      if (reload)
        [myrecording, opts] = track_spots(myrecording, opts);
        if (handles.autosave)
          save([myrecording.experiment '.mat'], 'myrecording', 'opts');
        end
      end
      set(handles.hFig, 'Visible', 'on')

    % Call the path filtering GUI, and filter accordingly
    case 'paths'
      set(handles.hFig, 'Visible', 'off')
      [myrecording, opts, reload] = inspect_paths(myrecording, opts);
      if (reload)
        [myrecording, opts] = filter_paths(myrecording, opts);
        if (handles.autosave)
          save([myrecording.experiment '.mat'], 'myrecording', 'opts');
        end
      end
      set(handles.hFig, 'Visible', 'on')
  end

  % Release the GUI and recompute the filters
  data.handles = handles;
  data.recording = myrecording;
  data.options = opts;
  guidata(handles.hFig, data);
  %if (reload)
  %  setup_environment()
  %end
  update_display(data, reload);
  set(handles.all_buttons, 'Enable', 'on');

  return
end

function experiment_Callback(hObject, eventdata)
% This function is responsible for handling the content of the
% structure which contains the parameters of the filtering algorithms.

  data = guidata(hObject);
  handles = data.handles;

  % Block the GUI
  %all_all = [handles.all_buttons, handles.save, handles.pipeline];
  %curr_status = get(all_all, 'Enable');
  %set(all_all, 'Enable', 'off');
  set(handles.all_buttons, 'Enable', 'off');
  drawnow;
  refresh(handles.hFig);

  % And get the type of button which called the callback (from its tag)
  type = get(hObject, 'tag');

  % By default, recompute
  recompute = true;
  nchannels = length(data.recording.channels);

  % Handle all three buttons differently
  switch type

    % Call the loading function
    case 'load'

      if (nchannels > 0)
        answer = questdlg('Save the current project ?', 'Save ?');
        if (strcmp(answer, 'Yes'))
          uisave({'myrecording','opts'}, [myrecording.experiment '.mat'])
        end
      end

      % Reset the current display
      if (all(ishandle(handles.img)))
        delete(handles.img);
        delete(handles.data);
        delete(handles.scale);
        dragzoom2D(handles.axes(1), 'off')
        set(handles.hFig, 'UserData', '');
      end

      % Fancy output
      disp('[Select a MAT file]');

      % Prompting the user for the MAT file
      [fname, dirpath] = uigetfile({'*.mat'}, ['Load a MAT file']);

      % Not cancelled
      if (ischar(fname))

        fname = fullfile(dirpath, fname);

        % Load the matrix and check its content
        matdata = load(fname);

        % Not what we expected
        if (~isfield(matdata, 'myrecording') || ~isfield(matdata, 'opts'))
          disp(['Error: ' fname ' does not contain a valid myrecording structure']);

        % Extract the loaded data
        else
          myrecording = matdata.myrecording;
          opts = update_structure(matdata.opts, 'options');

          data.recording = myrecording;
          data.options = opts;

          data.img = [];
          data.orig_img = [];
          data.img_next = [];
          data.handles.img = -1;
          data.handles.scatter = -1;
          data.handles.roi = {{}};
          data.handles.display = [1 1];
          data.handles.prev_channel = -1;
          data.handles.prev_frame = [-1 -1];
          data.handles.frame = [1 1];
          data.handles.prev_channel = -1;
          data.handles.current = 1;
        end
      end

    % Call the saving function
    case 'save'
      uisave({'myrecording','opts'}, [myrecording.experiment '.mat'])
      recompute = false;

    % Export the results
    case 'export'
      props = get_struct('exporting');
      props.file_name = myrecording.experiment;
      myrecording.channels = channels;

      [props, do_export] = edit_options(props);
      if (do_export)
        if (props.export_data)
          export_tracking(myrecording, props, opts);
        end
        if (props.export_movie)
          export_movie(myrecording, props, opts);
        end
      end
      recompute = false;
  end

  % Release the GUI and recompute the filters
  %set(all_all, {'Enable'}, curr_status);
  if (recompute)
  %  setup_environment()
  %  update_display(recompute);
    data = update_display(data, recompute);
    guidata(handles.hFig, data);
  end
  set(handles.all_buttons, 'Enable', 'on');

  return
end

function options_Callback(hObject, eventdata)
% This function is responsible for handling the buttons responsible for the
% option structure

  data = guidata(hObject);
  handles = data.handles;

  % Block the GUI
  set(handles.all_buttons, 'Enable', 'off');
  drawnow;
  refresh(handles.hFig);

  % And get the type of button which called the callback (from its tag)
  type = get(hObject, 'tag');

  % By default, recompute
  recompute = true;

  % Handle all three buttons differently
  switch type

    % Display the help text
    case 'help'
      % Dragzoom help message
      imghelp = regexp(help('dragzoom2D'), ...
                 '([ ]+Interactions:.*\S)\s+Adapted from','tokens');
      polyhelp = regexp(help('impoly'), ...
                 '([ ]+Interactions.*\S)\s+Simon','tokens');

      txthelp = ['AXES interactions:\n\n', imghelp{1}{1}, '\n\nROI interactions:\n\n', polyhelp{1}{1}];

      helpdlg(sprintf(txthelp), 'How to handle ASSET?');

    % Save a snapshot
    case 'snapshot'

      % Fancy output
      disp('[Select a SVG filename]');

      curdir = '';
      if(exist('export', 'dir'))
        curdir = pwd;
        cd('export');
      end

      % Prompting the user for the filename
      [fname, dirpath] = uiputfile({'*.svg', 'SVG vectorized image'}, ['Select a filename for your snapshot'], fullfile(pwd, 'snapshot.svg'));

      % Return back to our original folder
      if(~isempty(curdir))
        cd(curdir);
      end

      % Not cancelled
      if (ischar(fname))

        % This might take a while
        curr_name = get(handles.hFig, 'Name');
        set(handles.hFig, 'Name', [curr_name ' (Saving snapshot...)']);

        % Get the full name and save the snapshot !
        fname = fullfile(dirpath, fname);
        plot2svg(fname, handles.hFig);
        %print(handles.hFig, fname, '-dsvg');

        % And release !
        set(handles.hFig, 'Name', curr_name);
      end

      recompute = false;
  end

  % Release the GUI and recompute the display
  data = update_display(data, recompute);
  guidata(handles.hFig, data);
  set(handles.all_buttons, 'Enable', 'on');

  return
end

function gui_Callback(hObject, eventdata)
% This function handles the callback of most buttons in the GUI !

  data = guidata(hObject);
  handles = data.handles;

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
      data.recording.channels(indx).(type) = logical(get(hObject, 'Value'));

    % A change in the channel index
    case 'channels'
      handles.current = get(hObject, 'Value');

    % Otherwise, do nothing. This is used to cancel the deletion requests
    otherwise
      return;
  end

  data.handles = handles;

  % Update the display accordingly
  data = update_display(data, recompute);

  guidata(handles.hFig, data);

  return
end

%function cancel_CloseRequestFcn(hObject, eventdata)
%% This function stops the current processing after confirmation
%
%  data = guidata(hObject);
%  handles = data.handles;
%
%  % Just double check that the user want to quit
%  answer = questdlg('Do you really want to discard all your changes ?');
%  ok = strcmp(answer,'Yes');
%
%  % If everything is OK, release the GUI and quit
%  if (ok)
%    delete(handles.hFig);
%  end
%
%  return
%end

function gui_CloseRequestFcn(hObject, eventdata)
% This function converts the various indexes back into strings to prepare
% the channels structure for its standard form before releasing the GUI
% to exit it

  % If everything is OK, release the GUI and quit
  data = guidata(hObject);
  handles = data.handles;
  myrecording = data.recording;
  opts = data.options;

  if (strcmp(get(handles.hFig, 'Visible'), 'off'))
    ok = false;
    no = true;
  else
    % Just double check that the user want to save?
    answer = questdlg('Do you want to save your changes ?');
    ok = strcmp(answer,'Yes');
    no = strcmp(answer,'No');
  end

  % If everything is OK, release the GUI and quit
  if (ok)
    data = handle_rois(data, handles.current, 0);
    handles = data.handles;

    %if (handles.current > 0)
    %  data.recording.channels(handles.current).pixel_size = str2double(get(handles.resolution, 'String'));
    %  data.recording.channels(handles.current).amplitude = str2double(get(handles.amplitude, 'String'));
    %end

    % Copy the channels
    %myrecording.channels = data.recording.channels;
    % And get the experiment name
    myrecording.experiment = get(handles.experiment, 'String');

    [fpath, fname, fext] = fileparts(myrecording.channels(1).fname);
    mname = fullfile(fpath, [myrecording.experiment '.mat']);
    save(mname, 'myrecording');

    msgbox({'The analysis was saved to:', mname});

    delete(handles.hFig);
  elseif (no)
    delete(handles.hFig);
  end

  return
end

function fnames = load_images()

  % In case a subfolder name Movies exists, move into it for prompting
  curdir = '';
  if(exist('Movies', 'dir'))
    curdir = pwd;
    cd('Movies');
  elseif(exist(['..' filesep 'Movies'], 'dir'))
    curdir = pwd;
    cd(['..' filesep 'Movies']);
  end

  % Fancy output
  disp('[Select an embryo recording]');

  % Prompting the user for the movie file
  [fnames, dirpath] = uigetfile({'*.*'}, ['Load recordings of embryos, or one previously saved .MAT file'], 'MultiSelect', 'on');

  % Return back to our original folder
  if(~isempty(curdir))
    cd(curdir);
  end

  % If no file was selected, stop here
  if (isempty(fnames)  ||  isequal(dirpath, 0))
    disp(['No image selected']);
    fnames = '';
    return;
  end
  fnames = fullfile(dirpath, fnames);

  %% All to cells
  if (ischar(fnames))
    fnames = {fnames};
  end

  return;
end
