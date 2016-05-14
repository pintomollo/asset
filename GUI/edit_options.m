function [mystruct, is_updated] = edit_options(mystruct, name)
% EDIT_OPTIONS displays a GUI enabling the user to interactively
% modify the content of a structure provinding tooltip help from
% the comments of the corresponding section of get_struct.m.
%
%   [MODIFIED] = EDIT_OPTIONS(ORIGINAL) displays the content of
%   ORIGINAL to allow modifying it. Pressing the "Cancel" button
%   will ignore changes made by the user.
%
%   [...] = EDIT_OPTIONS(ORIGINAL, NAME) displays NAME in the title
%   bar. Used primarly during esition of sub-structures.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 09.05.14

  % Input check
  if (nargin == 1)
    name = '';
  end

  % Analyze the structure we received
  [new_name, values] = parse_struct(mystruct);
  if (isempty(name))
    name = new_name;
  end

  % Create the figure corresponding to the structure
  hFig = create_figure(values);
  is_updated = true;

  % Wait for the user to finish and delete the figure
  uiwait(hFig);
  delete(hFig);

  % Recompute the pixel size just in case
  mystruct = set_pixel_size(mystruct);

  return;

  % Function which creates the figure along with all the fields and controls
  function hFig = create_figure(myvals)

    % Fancy naming
    if (isempty(name))
      ftitle = 'Edit options';
    else
      ftitle = ['Edit options (' name ')'];
    end

    % THE figure
    hFig = figure('PaperUnits', 'centimeters',  ...
                  'CloseRequestFcn', @empty, ...            % Cannot be closed
                  'Color',  [0.7 0.7 0.7], ...
                  'MenuBar', 'none',  ...                   % No menu
                  'Name', ftitle,  ...
                  'Resize', 'off', ...                      % Cannot resize
                  'NumberTitle', 'off',  ...
                  'Units', 'normalized', ...
                  'Position', [0.3 0.25 .35 0.5], ...       % Fixed size
                  'DeleteFcn', @empty, ...                  % Cannot close
                  'WindowScrollWheelFcn', @scrolling, ...   % Supports scrolling
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

    % To accept changes
    hOK = uicontrol('Parent', hFig, ...
                  'Units', 'normalized',  ...
                  'Callback', @save_CloseRequestFcn, ...
                  'Position', [0.025 0.85 0.1 0.1], ...
                  'String', 'OK',  ...
                  'Tag', 'okbutton');

    % To discard changes
    hCancel = uicontrol('Parent', hFig, ...
                  'Units', 'normalized',  ...
                  'Callback', @cancel_CloseRequestFcn, ...
                  'Position', [0.025 0.75 0.1 0.1], ...
                  'String', 'Cancel',  ...
                  'Tag', 'okbutton');

    % Now we work in pixels, easier that way
    set(hPanel, 'Units', 'Pixels');

    % Get the original size of the panel, important to know the size
    % of the visible area.
    psize = get(hPanel, 'Position');

    % Storing all the handles for the controls, necessary to retrieve
    % their content
    fields = NaN(size(myvals, 1), 1);

    % We cycle through the list of fields and create the appropriate controls.
    % We start from the end of the structure to display it in the correct order.
    count = 0;
    for i=size(myvals,1):-1:1

      % Here is the trick, if we have drawn outside of the panel, increase and
      % and slide it !
      curr_size = get(hPanel, 'Position');
      if (count*50 + 70 > curr_size(4))
        set(hPanel, 'Position', curr_size+[0 -50 0 50]);
      end

      % The text defining the content of the field
      hText = uicontrol('Parent', hPanel, ...
                        'Units', 'pixels',  ...
                        'Position', [20 count*50 + 20 120 30], ...
                        'String', myvals{i,1},  ...
                        'Style', 'text',  ...
                        'TooltipString', myvals{i,5}, ...
                        'Tag', 'text');

      % The different types of controls used, this decision was made
      % in parse_struct.
      switch myvals{i,4}
        % Simple one-line editable fields
        case 'edit'
          hControl = uicontrol('Parent', hPanel, ...
                        'Units', 'pixels',  ...
                        'Position', [180 count*50 + 30 180 40], ...
                        'BackgroundColor', [1 1 1], ...
                        'String', myvals{i,2}, ...
                        'Style', myvals{i,4}, ...
                        'Tag', 'data');

        % Checkbox used for boolean values
        case 'checkbox'
          hControl = uicontrol('Parent', hPanel, ...
                        'Units', 'pixels',  ...
                        'Position', [180 count*50 + 30 180 40], ...
                        'Value', myvals{i,2}, ...
                        'Style', myvals{i,4}, ...
                        'Tag', 'data');

        % An incredibly flexible table, used for cell arrays
        case 'table'
          [nrows, ncolumns] = size(myvals{i,2});
          curr_size = [180 count*50 + 30 min(60*ncolumns + 5 + 30*(nrows>3), 180) min(20*nrows + 5 + 30*(ncolumns>3), 80)];
          hControl = uitable('Parent', hPanel, ...
                        'Units', 'pixels',  ...
                        'Position', curr_size, ...
                        'ColumnWidth', repmat({58}, 1, ncolumns), ...
                        'ColumnEditable', true(1, ncolumns), ...
                        'Data', myvals{i,2}, ...
                        'ColumnName', [], ...
                        'RowName', [], ...
                        'Tag', 'data');

          % Because of the sliders inherent to the table, it is too wide
          % so we increase its size and move the text a bit
          set(hText, 'Position', [20 count*50 + 10 + (curr_size(end)/2) 120 30]);
          %set(hText, 'Position', [20 count*50 + 50 120 30]);
          count = count + (curr_size(end)>50);

          if (strncmp(myvals{i,3}, 'strel', 5))
            set(hControl, 'Position', [180 (count-1)*50 + 30 180 120], ...
                          'ColumnWidth', repmat({20}, 1, ncolumns));
            set(hText, 'Position', [20 count*50 + 50 120 30]);

            count = count + 1;
          end

        % We use a simple text-field for ND cell arrays that cannot be edited
        case 'lock'
          hControl = uicontrol('Parent', hPanel, ...
                        'Units', 'pixels',  ...
                        'Position', [180 count*50 + 20 180 40], ...
                        'String', 'Sorry, cell arrays of cells cannot be edited',  ...
                        'Style', 'text',  ...
                        'Tag', myvals{i,1});

        % We use a push button to edit the underlying structure recursively
        case 'button'
          hControl = uicontrol('Parent', hPanel, ...
                        'Units', 'pixels',  ...
                        'Position', [180 count*50 + 30 180 40], ...
                        'String', 'edit structure', ...
                        'Callback', @recursive_edit, ...
                        'Style', 'pushbutton', ...
                        'Tag', myvals{i,1});
      end

      % We need to count how many items we display, and store their handlers
      count = count + 1;
      fields(i) = hControl;
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
                     'controls', fields, ...
                     'fix_offset', psize(2), ...
                     'slider', hSlider);

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

  % Update the panel according to the movements of the mouse scroll button
  function scrolling(hObject, scroll_struct)

    % Move a bit faster
    move = 2*scroll_struct.VerticalScrollCount * scroll_struct.VerticalScrollAmount;

    % Update the slider value
    handles = get(hFig, 'UserData');
    offsets = get(handles.slider, {'Value', 'Min', 'Max'});
    new_offset = offsets{1} - move;
    new_offset = max(new_offset, offsets{2});
    new_offset = min(new_offset, offsets{3});

    set(handles.slider,'Value', new_offset);

    % Update the panel position
    p = get(handles.panel, 'Position'); 
    set(handles.panel, 'Position',[p(1) -new_offset+handles.fix_offset p(3) p(4)])

    return;
  end

  % Recursively edit a structure
  function recursive_edit(hObject, eventdata, handles)

    % Hide the current panel
    set(hFig, 'Visible', 'off');
    fieldname = get(hObject, 'Tag');

    % Fancy naming for the panel title
    if (isempty(name))
      fname = fieldname;
    else
      fname = [name '.' fieldname];
    end

    % Edit the substructure
    [value, is_updated] = edit_options(mystruct.(fieldname), fname);

    % Store the new values
    mystruct.(fieldname) = value;

    % Display the old panel
    set(hFig, 'Visible', 'on');

    return
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

    handles = get(hFig, 'UserData');

    % Retrieve the different types of values that are displayed
    for i=1:size(values, 1)
      switch values{i,4}
        case 'edit'
          val = get(handles.controls(i), 'String');
        case 'checkbox'
          val = logical(get(handles.controls(i), 'Value'));
        case 'table'
          val = get(handles.controls(i), 'Data');
        otherwise
          continue;
      end

      % And convert them back to their original types
      if (~isempty(val))
        switch values{i,3}
          case 'cell'
            goods = ~cellfun('isempty', val);
            gcol = any(goods, 1);
            grow = any(goods, 2);
            val = val(grow, gcol);
          case 'num'
            [tmp, correct] = mystr2double(val);

            % Enforce the data type
            while (~correct)
              answer = inputdlg(['''' values{i,1} ''' is not a valid number, do you want to correct it ?'], 'Correct a numerical value', 1, {val});
              if (isempty(answer))
                [tmp, correct] = mystr2double(values{i, 2});
              else
                [tmp, correct] = mystr2double(answer{1});
              end
            end

            val = tmp;
          case 'bool'
            if (ischar(val))
              [tmp, correct] = mystr2double(val);

              % Enforce the data type
              while (~correct)
                answer = inputdlg(['''' values{i,1} ''' is not a valid number, do you want to correct it ?'], 'Correct a numerical value', 1, {val});
                if (isempty(answer))
                  [tmp, correct] = mystr2double(values{i, 2});
                else
                  [tmp, correct] = mystr2double(answer{1});
                end
              end

              val = logical(tmp);
            end
          case 'func'
            [tmp, correct] = mystr2func(val);

            % Enforce the data type
            while (~correct)
              if (iscell(val))
                val = char(val(~cellfun('isempty', val)));
                val = [val repmat(' ', size(val, 1), 1)];
                val = val.';
                val = strtrim(val(:).');
              end

              answer = inputdlg(['''' values{i,1} ''' is not a valid function, do you want to correct it ?'], 'Correct a function handle', 1, {val});
              if (isempty(answer))
                [tmp, correct] = mystr2func(values{i, 2});
              else
                [tmp, correct] = mystr2func(answer{1});
              end
            end

            val = tmp;

            if (length(val)==1)
              val = val{1};
            end
          case 'strel'
            val = strel('arbitrary', val);
        end
      else
        % The empty values keeping the proper type
        switch values{i,3}
          case 'cell'
            val = {{}};
          case 'num'
            val = NaN;
          case 'func'
            val = @(varargin)([]);
          case 'strel'
            val = strel('arbitrary', []);
        end
      end

      mystruct.(values{i,1}) = val;
    end

    % And resume
    uiresume(gcbf)

    return
  end
end

% Convert a string to a function, making some verifications in between
function [values, correct] = mystr2func(value)

  if (iscell(value))
    splits = value;
    splits = splits(~cellfun('isempty', splits));
  else
    splits = regexp(value, '\s+', 'split');
  end

  values = cellfun(@str2func, splits, 'UniformOutput', false);

  implicit = cellfun(@(x)(x(1) == '@'), splits);
  explicit = cellfun(@(x)(~isempty(which(x))), splits);
  correct = all(implicit | explicit);

  return;
end

% Convert a string to a list of double numbers
function [values, correct] = mystr2double(value)

  if (iscell(value))
    values = NaN(size(value));

    correct = true;
    for i=1:numel(value)
      [values(i), tmp] = mystr2double(value{i});
      correct = correct && tmp;
    end
  elseif (isempty(value))
    values = NaN;
    correct = true;
  elseif (ischar(value))
    splits = regexp(value, '\s+', 'split');
    values = str2double(splits);
    nans = cellfun(@(x)(strncmpi(x, 'nan', 3)), splits);
    correct = ~any(isnan(values) & ~nans);
  else
    values = value;
    correct = true;
  end

  return;
end

% Create the main table to store the values, their type and the type of uibutton to
% be displayed
function [name, values] = parse_struct(mystruct)

  % We get all fields
  fields = fieldnames(mystruct);
  values = cell(length(fields), 5);

  % And parse them
  for i=1:length(fields)
    field = fields{i};
    val = mystruct.(field);

    % First store the name and its value
    values{i, 1} = field;
    values{i, 2} = val;

    % Now check its type and act accordingly
    switch class(val)

      case {'double', 'single', 'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32', 'uint64'}
        values{i, 3} = 'num';
        if (isempty(val))
          values{i, 4} = 'edit';
          values{i, 2} = {''};
        else
          values{i, 4} = 'table';
          values{i, 2} = values{i,2};
        end
      case 'char'
        values{i, 3} = 'char';
        values{i, 4} = 'edit';
      case 'cell'
        if (~isempty(val) & strncmp(class(val{1}), 'function_handle', 15))
          values{i, 3} = 'func';
          values{i, 4} = 'table';

          tmp_cell = repmat({''}, size(val)+5);
          tmp_cell(1:size(val,1), 1:size(val,2)) = cellfun(@(x){func2str(x)}, val);
          values{i, 2} = tmp_cell;
        elseif (~any(cellfun('isclass', val, 'cell')))
          values{i, 3} = 'cell';
          values{i, 4} = 'table';

          tmp_cell = repmat({''}, size(val)+5);
          tmp_cell(1:size(val,1), 1:size(val,2)) = val;
          values{i, 2} = tmp_cell;
        else
          values{i, 3} = 'cell';
          values{i, 4} = 'lock';
        end
      case 'struct'
        values{i, 3} = 'struct';
        values{i, 4} = 'button';
      case 'logical'
        if (numel(values{i,2})>1)
          values{i, 3} = 'bool';
          values{i, 4} = 'table';
          values{i, 2} = values{i,2};
        else
          values{i, 3} = 'bool';
          values{i, 4} = 'checkbox';
          values{i, 2} = double(values{i, 2});
        end
      case 'function_handle'
        values{i, 3} = 'func';
        values{i, 4} = 'edit';
        values{i, 2} = func2str(values{i, 2});
      case 'strel'
        values{i, 3} = 'strel';
        values{i, 4} = 'table';
        values{i, 2} = getnhood(values{i, 2});
      otherwise
        values{i, 3} = class(val);
    end
  end

  % Get the help messages
  [name, helps] = extract_help_message(mystruct);
  [goods, indxs] = ismember(values(:,1), helps(:,1));
  values(goods, 5) = helps(indxs, 2);

  return;
end

function [name, curr_fields] = extract_help_message(mystruct)
% EXTRACT_HELP_MESSAGE reads the help message corresponding to the provided structure.
% The help messages are the corresponding comments at the end of the line defining
% the structure field in get_struct.m.
%
%   [NAME, HELP] = EXTRACT_HELP_MESSAGE(MYSTRUCT) returns the NAME of MYSTRUCT as found
%   in get_struct.m, and the HELP for all the fields found as a cell array.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 28.08.2014

  % Find the correct file
  fname = which('get_struct.m');

  % We open it in text mode 
  fid = fopen(fname,'rt');

  % Initialize the output values
  name = '';
  tmp_fields = cell(0, 2);
  curr_fields = tmp_fields;

  % The current fields in the structure we received, for comparison
  fields = fieldnames(mystruct);

  % We loop throughout the file, line by line
  line = fgetl(fid);
  skipping = true;
  while ischar(line)

    % We remove unsignificant white spaces
    line = strtrim(line);

    % We ignore empty lines and comment lines
    if (length(line) > 0 && line(1) ~= '%')

      % We first look for new structures
      tokens = regexp(line, 'case ''(.+)''','tokens');

      % We found one !
      if (~isempty(tokens))
        % If we have something stored, it means we found our structure, so stop here
        if (~isempty(curr_fields) & ~skipping)
          break;
        end

        % Otherwise, we store the current name and prepare the new values
        tmp_name = tokens{1}{1};
        curr_fields = tmp_fields;
        skipping = false;

      % If we are not skipping this whole structure, we can read the comments
      elseif (~skipping)

        % Get the field name and comments
        tokens = regexp(line, '[^'']*''([^,'']+)''\s*,[^%]*%\s*(.+)','tokens');

        % If we got something, we check whether it is part of the current structure
        if (~isempty(tokens))
          curr_field = tokens{1};

          % If not, we skip the entire structure
          if (~ismember(fields, curr_field{1}))
            skipping = true;

          % Otherwise we store it !
          else
            name = tmp_name;
            curr_fields(end+1,:) = tokens{1};
          end
        end
      end
    end

    % Process the next line
    line = fgetl(fid);
  end

  % And close the file
  fclose(fid);

  % Add empty help for the missing fields
  missing = ~ismember(fields, curr_fields(:,1));
  if (any(missing))
    curr_fields = [curr_fields; [fields(missing) cellstr(repmat(' ',sum(missing), 1))]];
  end

  return;
end
