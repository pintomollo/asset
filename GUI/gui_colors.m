function [mycolormap, is_updated] = gui_colors(indx)
% GUI_COLORS displays a GUI enabling the user to choose his favorite colormap.
%
%   [COLORMAP] = GUI_COLORS displays the GUI and returns the selected COLORMAP.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 28.08.14

  % Default value
  if (nargin == 0)
    indx = 1;
  end

  % Get the available colors, start with the first one
  colors = get_struct('colors');
  if (indx > length(colors.colormaps))
    indx = 1;
  end
  mycolormap = colors.colormaps{indx};

  % Create the figure corresponding to the structure
  hFig = create_figure(colors.colormaps);
  is_updated = true;

  % Wait for the user to finish and delete the figure
  uiwait(hFig);
  delete(hFig);

  % Find the index corresponding to the colormap
  if is_updated
    for i=1:length(colors.colormaps)
      if (strcmp(mycolormap, colors.colormaps{i}))
        mycolormap = i;
        break;
      end
    end
  else
    mycolormap = indx;
  end

  return;

  % Function which creates the figure along with all the fields and controls
  function hFig = create_figure(list)

    % Fancy naming
    ftitle = 'Select a color map';

    % THE figure
    hFig = figure('PaperUnits', 'centimeters',  ...
                  'CloseRequestFcn', @empty, ...            % Cannot be closed
                  'Color',  [0.7 0.7 0.7], ...
                  'MenuBar', 'none',  ...                   % No menu
                  'Name', ftitle,  ...
                  'Resize', 'off', ...                      % Cannot resize
                  'NumberTitle', 'off',  ...
                  'Units', 'normalized', ...
                  'Position', [0.3 0.25 0.35 0.3], ...       % Fixed size
                  'DeleteFcn', @empty, ...                  % Cannot close
                  'HandleVisibility', 'callback',  ...
                  'Tag', 'main_fig',  ...
                  'UserData', [], ...
                  'Visible', 'off');                        % Starts hidden

    % A slider used when the structure is too big to fit in the figure
    hSlider = uicontrol('Parent', hFig, ...
                  'Units', 'normalized',  ...
                  'Callback', @slider_Callback, ...
                  'Position', [0.15 0 0.05 1], ...
                  'Value', 0, ...
                  'Max', 1, ...
                  'Min', 0, ...
                  'Style', 'slider', ...
                  'Tag', 'slider1');

    % The panel which contains the details of the structure. This is essential
    % to be able to slide the controls around. Based on an idea picked up
    % online
    hPanel = uipanel('Parent', hFig, ...
                     'Title', '',  ...
                     'Units', 'normalized', ...
                     'Tag', 'uipanel',  ...
                     'Clipping', 'on',  ...
                     'Position', [0.2 0 0.8 1]);

    hRadio = uibuttongroup('Parent', hPanel, ...
                         'Units', 'normalized',  ...
                         'SelectionChangedFcn', @radio_Callback, ...
                         'Position', [0.05 0.05 0.2 0.1], ...
                         'tag', 'radio');

    % To accept changes
    hOK = uicontrol('Parent', hFig, ...
                  'Units', 'normalized',  ...
                  'Callback', @save_CloseRequestFcn, ...
                  'Position', [0.025 0.75 0.1 0.2], ...
                  'String', 'OK',  ...
                  'Tag', 'okbutton');

    % To discard changes
    hCancel = uicontrol('Parent', hFig, ...
                  'Units', 'normalized',  ...
                  'Callback', @cancel_CloseRequestFcn, ...
                  'Position', [0.025 0.55 0.1 0.2], ...
                  'String', 'Cancel',  ...
                  'Tag', 'okbutton');

    % Now we work in pixels, easier that way
    set(hPanel, 'Units', 'Pixels');
    set(hRadio, 'Units', 'Pixels');

    % Get the original size of the panel, important to know the size
    % of the visible area.
    psize = get(hPanel, 'Position');

    % Initialize the colors
    ncolors = 64;
    color_count = 0;
    all_axes = [];

    % We cycle through the list of colors and create the appropriate controls.
    % We start from the end of the structure to display it in the correct order.
    count = 0;
    for i=length(list):-1:1

      % Here is the trick, if we have drawn outside of the panel, increase and
      % and slide it !
      curr_size = get(hPanel, 'Position');
      if (count*50 + 70 > curr_size(4))
        set(hPanel, 'Position', curr_size+[0 -50 0 50]);

        for j = 1:length(all_axes)
          set(all_axes(j), 'Position', get(all_axes(j), 'Position') + [0 -1 0 0]);
        end
      end

      % The button and the text defining the colormap
      hControl = uicontrol('Parent', hRadio, ...
                         'Units', 'pixels',  ...
                         'Position', [30 count*50 + 5 120 30], ...
                         'Style', 'radiobutton',  ...
                         'Value', (i==indx), ...
                         'String', list{i},  ...
                         'Tag', ['radio' num2str(i)]);

      % Draw the axes
      hAxes = axes('Parent', hPanel, ...
                      'Units', 'pixels',  ...
                      'Position', [230 count*50 + 30 180 20], ...
                      'Visible', 'off', ...
                      'Tag', 'data');
      all_axes(end+1) = hAxes;

      % Create the image itself,adding a shift to create the full colormap
      hImg = image([1:ncolors], 'Parent', hAxes);

      % Adapt the display
      set(hAxes,'Visible', 'off');
      set(hRadio, 'Position', [20 20 200 count*50 + 40]);

      % Count the total number of colors
      colormap(hAxes, list{i});

      % We need to count how many items we display
      count = count + 1;
    end

    % We might need to resize the figure one last time
    curr_size = get(hPanel, 'Position');
    if (count*50 + 30 > curr_size(4))
      diff_size = curr_size(4) - (count*50 + 30);
      set(hPanel, 'Position', curr_size+[0 diff_size 0 -diff_size]);
    end

    % Remove the slider if not needed
    curr_size = get(hPanel, 'Position');
    if (curr_size(4) - psize(4) < 1)
      set(hSlider, 'Visible', 'off');
    else
      set(hSlider, 'Max', curr_size(4) - psize(4), 'Value', curr_size(4) - psize(4));
    end

    % Store everything and display the panel
    handles = struct('panel', hPanel, ...
                     'fix_offset', psize(2), ...
                     'slider', hSlider);

    % Display the panel
    set(hFig, 'UserData', handles, ...
              'Visible', 'on');

    return;
  end

  % Used to prevent the figure from closing
  function empty(hObject, eventdata, handles)
    return;
  end

  % Upon cancel, we simply exit without applying any modifications
  function cancel_CloseRequestFcn(hObject, eventdata, handles)

    is_updated = false;
    uiresume(gcbf)

    return;
  end

  % This function handles the callback of the radio buttons
  function radio_Callback(hObject, eventdata)

    tmp_tag = get(eventdata.NewValue, 'tag');
    mycolormap = colors.colormaps{str2double(tmp_tag(6:end))};

    return;
  end

  % Move using the slider values
  function slider_Callback(hObject, eventdata, handles)

    handles = get(hFig, 'UserData');

    % The slider value
    offset = get(handles.slider,'Value');

    % Update the panel position
    p = get(handles.panel, 'Position');
    set(handles.panel, 'Position',[p(1) -offset+handles.fix_offset p(3) p(4)])

    return
  end

  % Store the new values
  function save_CloseRequestFcn(hObject, eventdata, handles)

    % Just resume
    uiresume(gcbf)

    return
  end
end
