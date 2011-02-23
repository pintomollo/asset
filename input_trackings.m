function trackings = input_trackings(trackings, name)

mygray = [0:255]' / 255;
mygray = [mygray mygray mygray];

hFig = figure(...
'Units','characters',...
'Color',[0.701960784313725 0.701960784313725 0.701960784313725],...
'Colormap',mygray, ....
'MenuBar','none',...
'Name',['Manual Trackings Importer : ' name],...
'NumberTitle','off',...
'Position',[34 306 90.5 32],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','trackings_fig',...
'UserData',[],...
'Visible','off',...
'CloseRequestFcn',@trackings_fig_CloseRequestFcn,...
'DeleteFcn',@empty);

hGroups = uicontrol(...
'Parent',hFig,...
'Units','characters',...
'Callback',@list_groups_Select,...
'FontSize',10,...
'Position',[5 8 41 20],...
'String','',...
'Style','listbox',...
'Value',1,...
'Tag','list_groups');

hFiles = uicontrol(...
'Parent',hFig,...
'Units','characters',...
'FontSize',10,...
'Position',[48 6 37 22],...
'String','',...
'Style','listbox',...
'Value',1,...
'Tag','list_files');

hAddGroup = uicontrol(...
'Parent',hFig,...
'Units','characters',...
'Callback',@add_group,...
'FontSize',10,...
'Position',[10 3.5 12 2],...
'String','Add Group',...
'Tag','add_group');

hDeleteGroup = uicontrol(...
'Parent',hFig,...
'Units','characters',...
'Callback',@delete_group,...
'FontSize',10,...
'Position',[28 3.5 12 2],...
'String','Delete Group',...
'Tag','delete_group');

hAddFile = uicontrol(...
'Parent',hFig,...
'Units','characters',...
'Callback',@add_file,...
'FontSize',10,...
'Position',[53 3.5 12 2],...
'String','Add File',...
'Tag','add_file');

hDeleteFile = uicontrol(...
'Parent',hFig,...
'Units','characters',...
'Callback',@delete_file,...
'FontSize',10,...
'Position',[71.5 3.5 12 2],...
'String','Delete File',...
'Tag','delete_file');

hFileExpr = uicontrol(...
'Parent',hFig,...
'Units','characters',...
'FontSize',10,...
'Position',[5 6 30 2],...
'String','',...
'Style','edit',...
'Tag','file_expr');

hOK = uicontrol(...
'Parent',hFig,...
'Units','characters',...
'Callback',@trackings_fig_CloseRequestFcn,...
'FontSize',10,...
'Position',[41.5 0 14 2],...
'String','OK',...
'Tag','ok');

hText = uicontrol(...
'Parent',hFig,...
'Units','characters',...
'FontSize',10,...
'Position',[11 29 22 1.5],...
'String','Segmentation groups',...
'Style','text',...
'Tag','text1');

hText = uicontrol(...
'Parent',hFig,...
'Units','characters',...
'FontSize',10,...
'Position',[57 28.5 18 2],...
'String','Manual segmentations',...
'Style','text',...
'Tag','text2');

hLoadGroup = uicontrol(...
'Parent',hFig,...
'Units','characters',...
'Callback',@load_this_group,...
'FontSize',10,...
'Position',[35 6 11 2],...
'String','Load Group',...
'Tag','load_group');

handles = struct('uigroups', hGroups, ...
                 'uifiles', hFiles, ...
                 'uiexpr',hFileExpr,...
                 'trackings',trackings);

set(hFig, 'UserData', handles);
set(hFig,'Visible','on');

update_display(hFig, 1);             

uiwait(hFig);

handles = get(hFig, 'UserData');
trackings = handles.trackings;

delete(hFig);
drawnow;

function empty(hObject, eventdata, handles)

function update_display(hfig, indx)

  handles = get(hfig,'UserData');
  trackings = handles.trackings;

  if (~isfield(trackings,'child'))
    ngroups = 0;
  else
    ngroups = length(trackings.child);
  end
  groups = '';
  for i=1:ngroups
    groups = [groups trackings.child(i).name];

    if (i<ngroups)
      groups = [groups '|'];
    end
  end

  set(handles.uigroups, 'String', groups);

  if (isempty(indx) || indx <= 0 || indx > ngroups)
    indx = ngroups;
  end

  set(handles.uigroups, 'Value', indx);

  if (ngroups > 0)
    set(handles.uiexpr, 'String', trackings.child(indx).expr);

    nfiles = length(trackings.child(indx).files);
    files = '';
    for i=1:nfiles
      files = [files trackings.child(indx).files(i).fname];

      if (i<nfiles)
        files = [files '|'];
      end
    end

    set(handles.uifiles, 'String', files);
    set(handles.uifiles, 'Value', nfiles);
  end

% --- Executes on selection change in list_groups.
function list_groups_Select(hObject, eventdata, handles)
% hObject    handle to list_groups (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns list_groups contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_groups

  hfig = gcbf;
  handles = get(hfig,'UserData');
  indx = get(handles.uigroups,'Value');

  update_display(hfig, indx);

% --- Executes on button press in add_group.
function add_group(hObject, eventdata, handles)
% hObject    handle to add_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  hfig = gcbf;
  handles = get(hfig,'UserData');
  trackings = handles.trackings;

  curdir = '';
  if(exist('Shapes', 'dir'))
    curdir = cd;
    cd('Shapes');
  elseif(exist('../Shapes', 'dir'))
    curdir = cd;
    cd('../Shapes');
  end

  [fname,dirpath] = uigetfile({'*.shapes', '*.mat'}, 'Load tracking(s) file', 'MultiSelect', 'on');

  if(~isempty(curdir))
    cd(curdir);
  end

  if (length(fname)==0  ||  isequal(dirpath,0))
    return;
  elseif (~iscell(fname))
    fname = {fname};
  end  

  for f=1:length(fname)
    [tokens,junk]=regexp(fname{f},'(.+[-_])?([^-_\.]+)(\..+)','tokens');
    name = tokens{1}{1};
    suffix = tokens{1}{2};
    ext = tokens{1}{3};

    if (strncmp(ext, '.mat', 4))
      tmp = load([dirpath fname{f}]);
      tmp_trackings = tmp.trackings;

      for i=1:length(tmp_trackings.child)
        found = false;
        for j=1:length(trackings.child)
          if (strcmp(trackings.child(j).expr, tmp_trackings.child(i).expr))
            found = true;
            break;
          end
        end

        if (~found)
          trackings.child(end+1) = tmp_trackings.child(i)
        end
      end
    else
      if (isempty(name))
        expr = [suffix ext];
        name = suffix;
      else
        expr = [name '*' ext];
      end

      found = false;
      for i=1:length(trackings.child)
        if (strcmp(trackings.child(i).expr, [dirpath expr]))
          found = true;
          break;
        end
      end

      if (~found)
        if (length(trackings.child) == 0)
          trackings.child = get_struct('tracking');
        end
        trackings.child(end+1).name = name;
        trackings.child(end).expr = [dirpath expr];
        trackings.child(end) = load_group(trackings.child(end));
      end
    end
  end

  handles.trackings = trackings;
  set(hfig,'UserData',handles)
  update_display(hfig, 0);

% --- Executes on button press in delete_group.
function delete_group(hObject, eventdata, handles)
% hObject    handle to delete_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hfig = gcbf;
handles = get(hfig, 'UserData');
trackings = handles.trackings;
indx = get(handles.uigroups,'Value');

if (indx > 0 && indx <= length(trackings.child))
  trackings.child(indx) = [];

  handles.trackings = trackings;
  set(hfig, 'UserData', handles);
  update_display(hfig, indx);
end

% --- Executes on button press in add_file.
function add_file(hObject, eventdata, handles)
% hObject    handle to add_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  hfig = gcbf;
  handles = get(hfig,'UserData');
  trackings = handles.trackings;
  indx = get(handles.uigroups,'Value');

  if (indx <= 0 || indx > length(trackings.child))
    hwarn = warndlg('You need to select a tracking group first', 'Add file Error');
    %setfocus(hwarn);

    return;
  end

  curdir = '';
  if(exist('Shapes', 'dir'))
    curdir = cd;
    cd('Shapes');
  elseif(exist('../Shapes', 'dir'))
    curdir = cd;
    cd('../Shapes');
  end

  [fname,dirpath] = uigetfile({'*.shapes'}, 'Load tracking(s) file', 'MultiSelect', 'on');

  if(~isempty(curdir))
    cd(curdir);
  end

  if (length(fname)==0  ||  isequal(dirpath,0))
    return;
  elseif (~iscell(fname))
    fname = {fname};
  end  

  for f=1:length(fname)
    found = false;
    for k=1:length(trackings.child(indx).files)
      if (strcmp([dirpath fname{f}], trackings.child(indx).files(k).fname))
        found = true;
        break;
      end
    end
    if (~found)
      if (length(trackings.child(indx).files) == 0)
        trackings.child(indx).files = get_struct('file');
      end
      trackings.child(indx).files(end+1).fname = [dirpath fname{f}];
    end
  end

  handles.trackings = trackings;
  set(hfig,'UserData',handles)
  update_display(hfig, indx);

% --- Executes on button press in delete_file.
function delete_file(hObject, eventdata, handles)
% hObject    handle to delete_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  hfig = gcbf;
  handles = get(hfig,'UserData');
  trackings = handles.trackings;
  gindx = get(handles.uigroups,'Value');
  findx = get(handles.uifiles,'Value');

  if (gindx > 0 && findx > 0 && findx <= length(trackings.child(gindx).files))
    trackings.child(gindx).files(findx) = [];
    handles.trackings = trackings;
    set(hfig, 'UserData', handles);
  end

  update_display(hfig, gindx);

% --- Executes on button press in load_group.
function load_this_group(hObject, eventdata, handles)
% hObject    handle to load_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  hfig = gcbf;
  handles = get(hfig,'UserData');
  trackings = handles.trackings;
  indx = get(handles.uigroups,'Value');
  expr = get(handles.uiexpr, 'String');

  if (indx <= 0 || indx > length(trackings.child))
    return;
  end

  if (~strcmp(expr, trackings.child(indx).expr))
    trackings.child(indx).expr = expr;
  end

  trackings.child(indx) = load_group(trackings.child(indx));

  handles.trackings = trackings;
  set(hfig, 'UserData', handles);

  update_display(hfig, indx);

function trackings = load_group(trackings)

  expr = trackings.expr;
  
  slash = findstr(expr, '/');
  dirpath  = '';
  if (length(slash)>0)
    dirpath = expr(1:slash(end));
  end

  files = dir(expr);
  if (length(files) == 0)
    hwarn = warndlg('No file correspond to this regular expression', 'Load files Error');
    %setfocus(hwarn);
  else
    for f = 1:length(files)
      name = [dirpath files(f).name];
      found = false;
      for k=1:length(trackings.files)
        if (strcmp(name, trackings.files(k).fname))
          found = true;
          break;
        end
      end
      if (~found)
        if (length(trackings.files) == 0)
          trackings.files = get_struct('file');
        end
        trackings.files(end+1).fname = name;
      end
    end
  end

% --- Executes on button press in ok.
function trackings_fig_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hfig = gcbf;
uiresume(hfig);

