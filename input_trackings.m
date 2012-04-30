function trackings = input_trackings(trackings, name, opts)
% INPUT_TRACKING displays a pop-up window to import your .shapes tracking files.
%
%   TRACKINGS = INPUT_TRACKINGS(TRACKINGS, OPTS) fills TRACKINGS with the corresponding
%   data from the files (see load_shapes.m) that are provided interactively by the user.
%   
%   [...] = INPUT_TRACKING(..., NAME) displays NAME in the title bar as a reminder for
%   the current experiment which the trackings correspond to.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 20.03.2011

  % Some basic input checking to handle missing arguments
  if (nargin < 2)
    name = '';
    opts = get_struct('ASSET');
  elseif (nargin < 3)
    if (isstruct(name))
      opts = name;
      name = '';
    else
      opts = get_struct('ASSET');
    end
  end

  % If trackings does not have the correct structure, merge it with the standard one
  if (~isfield(trackings, 'child'))
    trackings = merge_structures(trackings, get_struct('tracking'));
  end

  % Build our own RGB grayscale colormap
  mygray = [0:255]' / 255;
  mygray = [mygray mygray mygray];

  % Create the main window, non-resizable, non-closeable 
  hFig = figure( ...
      'Units', 'characters', ...
      'Color', [0.7 0.7 0.7], ...
      'Colormap', mygray, ....
      'MenuBar', 'none', ...
      'Name', ['Manual Trackings Importer : ' name], ...
      'NumberTitle', 'off', ...
      'Position', [34 306 90.5 32], ...
      'Resize', 'off', ...
      'HandleVisibility', 'callback', ...
      'Tag', 'trackings_fig', ...
      'UserData', [], ...
      'Visible', 'off', ...
      'CloseRequestFcn', @close_trackings_fig, ...
      'DeleteFcn', @empty);

  % Create the file groups list box (shows all the select groups of files)
  hGroups = uicontrol( ...
      'Parent', hFig, ...
      'Units', 'characters', ...
      'Callback', @list_groups_select, ...
      'FontSize', 10, ...
      'Position', [5 8 41 20], ...
      'String', '', ...
      'Style', 'listbox', ...
      'Value', 1, ...
      'Tag', 'list_groups');

  % The text above the goup list box
  hText = uicontrol( ...
      'Parent', hFig, ...
      'Units', 'characters', ...
      'FontSize', 10, ...
      'Position', [11 29 22 1.5], ...
      'String', 'Segmentation groups', ...
      'Style', 'text', ...
      'Tag', 'text1');

  % Create the files list box (shows the individial files in the current group)
  hFiles = uicontrol( ...
      'Parent', hFig, ...
      'Units', 'characters', ...
      'FontSize', 10, ...
      'Position', [48 6 37 22], ...
      'String', '', ...
      'Style', 'listbox', ...
      'Value', 1, ...
      'Tag', 'list_files');

  % The text above the file list box
  hText = uicontrol( ...
      'Parent', hFig, ...
      'Units', 'characters', ...
      'FontSize', 10, ...
      'Position', [57 28.5 18 2], ...
      'String', 'Manual segmentations', ...
      'Style', 'text', ...
      'Tag', 'text2');

  % The button to add a new file group
  hAddGroup = uicontrol( ...
      'Parent', hFig, ...
      'Units', 'characters', ...
      'Callback', @add_group, ...
      'FontSize', 10, ...
      'Position', [10 3.5 12 2], ...
      'String', 'Add Group', ...
      'Tag', 'add_group');

  % The button used to load the group defined by the current regular expression
  hLoadGroup = uicontrol( ...
      'Parent', hFig, ...
      'Units', 'characters', ...
      'Callback', @load_this_group, ...
      'FontSize', 10, ...
      'Position', [35 6 11 2], ...
      'String', 'Load Group', ...
      'Tag', 'load_group');

  % The button to delete a file group
  hDeleteGroup = uicontrol( ...
      'Parent', hFig, ...
      'Units', 'characters', ...
      'Callback', @delete_group, ...
      'FontSize', 10, ...
      'Position', [28 3.5 12 2], ...
      'String', 'Delete Group', ...
      'Tag', 'delete_group');

  % The button to add an individual file to the current group
  hAddFile = uicontrol( ...
      'Parent', hFig, ...
      'Units', 'characters', ...
      'Callback', @add_file, ...
      'FontSize', 10, ...
      'Position', [53 3.5 12 2], ...
      'String', 'Add File', ...
      'Tag', 'add_file');

  % The button to delete an individual file from the current group
  hDeleteFile = uicontrol( ...
      'Parent', hFig, ...
      'Units', 'characters', ...
      'Callback', @delete_file, ...
      'FontSize', 10, ...
      'Position', [71.5 3.5 12 2], ...
      'String', 'Delete File', ...
      'Tag', 'delete_file');

  % The regular expression defining the current file group
  hFileExpr = uicontrol( ...
      'Parent', hFig, ...
      'Units', 'characters', ...
      'FontSize', 10, ...
      'Position', [5 6 30 2], ...
      'String', '', ...
      'Style', 'edit', ...
      'Tag', 'file_expr');

  % The OK button to quit the interface
  hOK = uicontrol( ...
      'Parent', hFig, ...
      'Units', 'characters', ...
      'Callback', @close_trackings_fig, ...
      'FontSize', 10, ...
      'Position', [41.5 0 14 2], ...
      'String', 'OK', ...
      'Tag', 'ok');

  % The structure stored in the UserData to pass arguments to the callback functions
  handles = struct('uigroups', hGroups, ...
                   'uifiles', hFiles, ...
                   'uiexpr', hFileExpr, ...
                   'opts', opts, ...
                   'trackings', trackings);

  % Stores the arguments in UserData and show the GUI
  set(hFig, 'UserData', handles);
  set(hFig, 'Visible', 'on');

  % Update the window
  update_display(hFig, 1);             

  % Wait for the window to quit
  uiwait(hFig);

  % Retrieve the updated trackings structure
  handles = get(hFig, 'UserData');
  trackings = handles.trackings;

  % Destroy the figure now
  delete(hFig);
  drawnow;

  return;

  % This function is used only to bypass the window closing so that we make sure
  % to return the trackings structure
  function empty(hObject, eventdata, handles)

    return;
  end

  % This function updates the display of the window, to include new groups and
  % show the files of the corresponding group
  function update_display(hfig, indx)

    % Get the data
    handles = get(hfig,'UserData');
    trackings = handles.trackings;

    % Count the number of groups currently in trackings
    if (~isfield(trackings,'child'))
      ngroups = 0;
    else
      ngroups = length(trackings.child);
    end

    % Prepare the string used to buil the group list box
    groups = '';

    % Loop over each group
    for i=1:ngroups

      % Concatenate its name to the growing list
      groups = [groups trackings.child(i).name];

      % Separate them with | such that it create distinct elements of the list
      if (i<ngroups)
        groups = [groups '|'];
      end
    end

    % Update the content of the group list
    set(handles.uigroups, 'String', groups);

    % If the index of the current group is not valid, choose the last one
    if (isempty(indx) || indx <= 0 || indx > ngroups)
      indx = ngroups;
    end

    % Select the corresponding element in the list
    set(handles.uigroups, 'Value', indx);

    % If we are actually selecting a group, build the corresponding list of files the
    % same way
    if (ngroups > 0)
      
      % Update also the regular expression
      set(handles.uiexpr, 'String', trackings.child(indx).expr);

      nfiles = length(trackings.child(indx).child);
      files = '';
      for i=1:nfiles
        files = [files trackings.child(indx).child(i).fname];

        if (i<nfiles)
          files = [files '|'];
        end
      end

      set(handles.uifiles, 'String', files);
      set(handles.uifiles, 'Value', nfiles);
    end

    return;
  end

  % Retrieves the selected group and update the display
  function list_groups_select(hObject, eventdata, handles)

    hfig = gcbf;
    handles = get(hfig,'UserData');
    indx = get(handles.uigroups,'Value');

    update_display(hfig, indx);

    return;
  end

  % Add a group of files to trackings
  function add_group(hObject, eventdata, handles)

    % Get the UserData and the trackings
    hfig = gcbf;
    handles = get(hfig,'UserData');
    trackings = handles.trackings;
    opts = handles.opts;

    % Start in the current directory
    curdir = '';

    % If a Shapes sub-directory exists, use it instead
    if(exist('Shapes', 'dir'))
      
      % Update the current directory accordingly and move into Shapes
      curdir = cd;
      cd('Shapes');

    % Otherwise, check in the parent directory
    elseif(exist(['..' filesep 'Shapes'], 'dir'))
      curdir = cd;
      cd(['..' filesep 'Shapes']);
    end

    % Open a pop-up to select .shapes file(s)
    [fname, dirpath] = uigetfile({'*.shapes', '*.mat'}, 'Load tracking(s) file', 'MultiSelect', 'on');

    % If we changed directory, go back to our initial one
    if(~isempty(curdir))
      cd(curdir);
    end

    % If no file was selected, exit
    if (length(fname) == 0  ||  isequal(dirpath, 0))
      return;
    
    % If only one file was selected, treat it as a single cell
    elseif (~iscell(fname))
      fname = {fname};
    end  

    % Loop over the selected files
    for f = 1:length(fname)

      % Extract the name/suffix/extension from the filename according to the provided regular epxression
      [tokens, junk] = regexp(fname{f}, opts.file_regexpr, 'tokens');
      name = tokens{1}{1};
      suffix = tokens{1}{2};
      ext = tokens{1}{3};

      % If a .mat file was selected, load it
      if (strncmp(ext, '.mat', 4))
        tmp = load([dirpath fname{f}]);

        % If it do not contain a valid trackings variable, skip this file
        if (~isefield(tmp, 'trackings') | ~isefield(tmp.trackings, 'child'))
          continue;
        end

        % Get the trackings variable
        tmp_trackings = tmp.trackings;

        % Loop over the children of the trackings
        for i=1:length(tmp_trackings.child)

          % For now on we have not found the same group in our list
          found = false;

          % Loop over our trackings
          for j=1:length(trackings.child)

            % If the regular expression is the same, we already have this group
            if (strcmp(trackings.child(j).expr, tmp_trackings.child(i).expr))
              found = true;
              break;
            end
          end

          % If this group does not exist, add it to the list
          if (~found)
            trackings.child(end+1) = tmp_trackings.child(i)
          end
        end

      % We have selected a normal file
      else

        % Create the regular expression
        if (isempty(name))
          expr = [suffix ext];
          name = suffix;
        else
          expr = [name '*' ext];
        end

        % Check if we already have it in our list
        found = false;
        for i=1:length(trackings.child)
          tmp_expr = [dirpath expr];
          if (strcmp(trackings.child(i).expr, tmp_expr) | strcmp(trackings.child(i).expr, relativepath(tmp_expr)))

            found = true;
            break;
          end
        end

        % Otherwise, add it to our list of groups
        if (~found)
          
          % Handle empty children
          if (length(trackings.child) == 0)
            trackings.child = get_struct('tracking', 0);
          end

          % Store the required data
          trackings.child(end+1).name = name;
          trackings.child(end).expr = [dirpath expr];
          trackings.child(end) = load_group(trackings.child(end));
        end
      end
    end

    % Store the data and update the display
    handles.trackings = trackings;
    set(hfig,'UserData',handles)
    update_display(hfig, 0);

    return;
  end

  % Deletes the selected group of files
  function delete_group(hObject, eventdata, handles)

    % Retrieve the necessary data
    hfig = gcbf;
    handles = get(hfig, 'UserData');
    trackings = handles.trackings;
    indx = get(handles.uigroups,'Value');

    % If the index is valid, delete the corresponding group
    if (indx > 0 && indx <= length(trackings.child))
      trackings.child(indx) = [];

      % Update the data and the display
      handles.trackings = trackings;
      set(hfig, 'UserData', handles);
      update_display(hfig, indx);
    end

    return;
  end

  % Add a file the current group
  function add_file(hObject, eventdata, handles)

    % Get the data
    hfig = gcbf;
    handles = get(hfig,'UserData');
    trackings = handles.trackings;
    indx = get(handles.uigroups,'Value');

    % If the index is not valid, explain it
    if (indx <= 0 || indx > length(trackings.child))
      warndlg('You need to select a tracking group first', 'Add file Error');

      return;
    end

    % Check whether we can find a Shapes dir and move into it if possible
    curdir = '';
    if(exist('Shapes', 'dir'))
      curdir = cd;
      cd('Shapes');
    elseif(exist(['..' filesep 'Shapes'], 'dir'))
      curdir = cd;
      cd(['..' filesep 'Shapes']);
    end

    % Open a pop-up to get the files
    [fname, dirpath] = uigetfile({'*.shapes'}, 'Load tracking(s) file', 'MultiSelect', 'on');

    % Move back to the original directory
    if(~isempty(curdir))
      cd(curdir);
    end

    % If no file were selected, exit
    if (length(fname)==0  ||  isequal(dirpath,0))
      return;

    % If only one file was selected, create a one cell array
    elseif (~iscell(fname))
      fname = {fname};
    end  

    % Loop over the files
    for f=1:length(fname)

      % Check whether we already have it in the current group
      found = false;
      for k=1:length(trackings.child(indx).child)

        % If we do, stop our search
        if (strcmp([dirpath fname{f}], trackings.child(indx).child(k).fname))
          found = true;

          break;
        end
      end

      % If we do not have it, add it to the current group
      if (~found)

        % Handle empty children
        if (length(trackings.child(indx).child) == 0)
          trackings.child(indx).child = get_struct('file', 0);
        end

        % Store the file
        trackings.child(indx).child(end+1).fname = [dirpath fname{f}];
      end
    end

    % Store the data and update the display
    handles.trackings = trackings;
    set(hfig,'UserData',handles)
    update_display(hfig, indx);

    return;
  end

  % Delete the current file
  function delete_file(hObject, eventdata, handles)

    % Get the data
    hfig = gcbf;
    handles = get(hfig,'UserData');
    trackings = handles.trackings;

    % Retrieve the index of group and file selected
    gindx = get(handles.uigroups,'Value');
    findx = get(handles.uifiles,'Value');

    % If both indexes are valid, delete the selected files
    if (gindx > 0 && findx > 0 && findx <= length(trackings.child(gindx).child))
      trackings.child(gindx).child(findx) = [];

      % Store the data and update the display
      handles.trackings = trackings;
      set(hfig, 'UserData', handles);
      update_display(hfig, gindx);
    end

    return;
  end

  % Load the group defined by the regular expression
  function load_this_group(hObject, eventdata, handles)
  
    % Retrieve all the necessary data
    hfig = gcbf;
    handles = get(hfig,'UserData');
    trackings = handles.trackings;
    indx = get(handles.uigroups,'Value');
    expr = get(handles.uiexpr, 'String');

    % If there are currently no groups, create the first one
    if (length(trackings.child) == 0)
      trackings.child = get_struct('tracking', 0);
      tmp_tracking = get_struct('tracking');
      indx = 1;

    % If the index is not valid, exit
    elseif (indx <= 0 || indx > length(trackings.child))
    
      return;

    % Otherwise, retrieve the selected group
    else
      tmp_tracking = trackings.child(indx);
    end

    % Update the regular expression if required
    if (~strcmp(expr, tmp_tracking.expr))
      tmp_tracking.expr = expr;
    end

    % Load the corresponding group and update trackings
    trackings.child(indx) = load_group(tmp_tracking);

    % Store the data and update the display
    handles.trackings = trackings;
    set(hfig, 'UserData', handles);
    update_display(hfig, indx);

    return;
  end

  % Load a group based on its regular expression
  function trackings = load_group(trackings)

    % Get the regular expression
    expr = trackings.expr;
    if (length(expr) < 1)
      return;
    end

    % Create an absolute path to the file
    if (expr(1) ~= filesep)
      expr = fullfile(pwd, expr);
    end
    
    % Extract only the name of the file
    slash = findstr(expr, filesep);
    dirpath  = '';
    if (length(slash)>0)
      dirpath = expr(1:slash(end));
      expr_name = expr(slash(end)+1:end);
    end

    % List the files on the hard disk
    files = dir(expr);

    % If we have not found any, notify the user
    if (length(files) == 0)
      warndlg('No file correspond to this regular expression', 'Load files Error');

    % Otherwise, parse them
    else
      common_string = '';

      % Loop over the results
      for f = 1:length(files)
        % Create the absolute path
        name = [dirpath files(f).name];

        % Now check if we have it already or not in our trackings
        found = false;
        for k=1:length(trackings.child)
          if (strcmp(name, trackings.child(k).fname))
            found = true;
            break;
          end
        end

        % If we don't, add it to the current group of trackings
        if (~found)
          if (length(trackings.child) == 0)
            trackings.child = get_struct('file', 0);
          end
          trackings.child(end+1).fname = name;

          % Extract the common part of the file names to create our experiment name
          if (isempty(common_string))
            common_string = files(f).name;
          else
            common_string = common_substring(common_string, files(f).name);
          end
        end
      end

      % Store the name
      if (isempty(common_string))
        trackings.name = regexprep(expr_name, '^\W+|\W+$|\*', '');
      else
        trackings.name = regexprep(common_string, '^\W+|\W+$', '');
      end
    end

    return;
  end

  % Closes the figure and resumes the main function, returning the updated trackings 
  function close_trackings_fig(hObject, eventdata, handles)

    hfig = gcbf;
    uiresume(hfig);

    return;
  end
end
