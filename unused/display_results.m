function display_results(mymovie, trackings, nframe)

if (nargin == 1)
  trackings = [];
  nframe = 1;
elseif (nargin == 2)
  if (isnumeric(trackings))
    nframe = trackings;
    trackings = [];
  else
    nframe = 1;
  end
end

mygray = [0:255]' / 255;
mygray = [mygray mygray mygray];

hFig = figure(...
'Units','characters',...
'Color',[0.701960784313725 0.701960784313725 0.701960784313725],...
'Colormap',mygray, ....
'KeyPressFcn',@KeyPressFcn, ...
'MenuBar','none',...
'NumberTitle','off',...
'Position',[0 0 212 54],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','trackings_fig',...
'UserData',[],...
'Visible','off');

hEggshell = axes(...
'Parent',hFig,...
'Units','characters',...
'Position',[3.5 9.5 95 41],...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'FontSize',10,...
'LooseInset',[27.5 6 20 4],...
'Tag','eggshell');

hCortex = axes(...
'Parent',hFig,...
'Units','characters',...
'Position',[110 9.5 95 41],...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'FontSize',10,...
'LooseInset',[27.5 6 20 4],...
'Tag','cortex');

linkaxes([hEggshell,hCortex]);

hType = uicontrol(...
'Parent',hFig,...
'Units','characters',...
'Callback',@change_type,...
'KeyPressFcn',@KeyPressFcn, ...
'FontSize',10,...
'Position',[96.5 50 20 2.5],...
'String', 'dic|markers|data',...
'Style','popupmenu',...
'Value',1,...
'Tag','type');

hText = uicontrol(...
'Parent',hFig,...
'KeyPressFcn',@KeyPressFcn, ...
'Units','characters',...
'FontSize',10,...
'Position',[35.5 51 27 2.5],...
'String','Eggshell',...
'Style','text',...
'Tag','text1');

hText = uicontrol(...
'Parent',hFig,...
'KeyPressFcn',@KeyPressFcn, ...
'Units','characters',...
'FontSize',10,...
'Position',[144 51 27 2.5],...
'String','Cortex',...
'Style','text',...
'Tag','text2');

hPanel = uipanel(...
'Parent',hFig,...
'Units','characters',...
'FontSize',10,...
'Title','Errors',...
'Tag','uipanel',...
'Clipping','on',...
'Position',[4 1 205 8]);

hTable1 = uitable(...
'Parent',hPanel,...
'KeyPressFcn',@KeyPressFcn, ...
'Units','characters',...
'Data',{  },...
'Position',[1 0.5 90.5 6.5],...
'Tag','uieggshell');

hTable2 = uitable(...
'Parent',hPanel,...
'KeyPressFcn',@KeyPressFcn, ...
'Units','characters',...
'Data',{},...
'Position',[115 0.5 90.5 6.5],...
'Tag','uicortex');

handles = struct('trackings',trackings, ...
                 'mymovie',mymovie, ...
                 'eggshell',struct('axes', hEggshell, ...
                                   'img', -1, ...
                                   'lines', [], ...
                                   'pts',[], ...
                                   'ruffles', []), ...
                 'cortex',struct('axes', hCortex, ...
                                 'img', -1, ...
                                 'lines', [], ...
                                 'pts', [], ...
                                 'ruffles', []), ...
                 'current', nframe, ...
                 'previous', nframe, ...
                 'type', 'dic', ...
                 'uitype',hType);

%                                   'errors', hTable1), ...
%                                 'errors', hTable2), ...

set(hFig, 'UserData', handles);
update_display(hFig);

set(hFig,'Visible','on');
end

function KeyPressFcn(hobj,event_data)
% Callback for key press in figure window.

% Get the ASCII value key.

key = double(event_data.Character);

% Shift, Alt, and Ctrl keys, when held down, generate 
% callbacks at the keyboard repeat rate, and these keys
% by themselves are stored as ''.

if isempty(key)
  return;
end   

modifier = event_data.Modifier;
%if (length(modifier) > 0)
%&&  strcmp(modifier{1},'control')

%  modifier

%end

hfig = gcbf;
handles = get(hfig,'UserData');
if (isfield(handles.mymovie, 'dic') & ~isempty(handles.mymovie.dic))
  [nframes] = size_data(handles.mymovie.dic);
elseif (isfield(handles.mymovie, 'cortex') & ~isempty(handles.mymovie.cortex))
  [nframes] = size_data(handles.mymovie.cortex);
elseif (isfield(handles.mymovie, 'eggshell') & ~isempty(handles.mymovie.eggshell))
  [nframes] = size_data(handles.mymovie.eggshell);
elseif (isfield(handles.mymovie, 'data') & ~isempty(handles.mymovie.data))
  [nframes] = size_data(handles.mymovie.data);
end

moved = false;

switch key
  case 28

    if (isempty(modifier))
      handles.previous = handles.current;

      handles.current = handles.current - 1;
      if (handles.current < 1)
        handles.current = nframes;
      end

      moved = true;
    else

      xaxes = get(handles.eggshell.axes,'XLim');
      xaxes = xaxes - (diff(xaxes) /4);

      xaxes = xaxes - (xaxes(1) * (xaxes(1) < 0));

      set(handles.eggshell.axes,'XLim',xaxes);
    end

% Left Arrow keypress:  rotate camera clockwise (as viewed
% from positive z).

  case 29

% Right Arrow keypress:  rotate camera ccw.

  if (isempty(modifier))
      handles.previous = handles.current;
      handles.current = handles.current + 1;
      if (handles.current > nframes)
        handles.current = 1;
      end

      moved = true;
  else
      maxx = size(get(handles.eggshell.img,'CData'),2);

      xaxes = get(handles.eggshell.axes,'XLim');
      xaxes = xaxes + (diff(xaxes) /4);

      xaxes = (xaxes - (xaxes(2) - maxx) * (xaxes(2) > maxx)) ;

      set(handles.eggshell.axes,'XLim',xaxes);
    end

case 30

% Up Arrow keypress:  increase elevation.
  
  if (~isempty(modifier))

    yaxes = get(handles.eggshell.axes,'YLim');
    yaxes = yaxes - (diff(yaxes) /4);
      
    yaxes = yaxes - (yaxes(1) * (yaxes(1) < 0));

    set(handles.eggshell.axes,'YLim',yaxes);
  end

case 31

% Down Arrow keypress:  decrease elevation.
  if (~isempty(modifier))
    maxy = size(get(handles.eggshell.img,'CData'),1);

    yaxes = get(handles.eggshell.axes,'YLim');
    yaxes = yaxes + (diff(yaxes) /4);

    yaxes = (yaxes - (yaxes(2) - maxy) * (yaxes(2) > maxy));

    set(handles.eggshell.axes,'YLim',yaxes);
  end

case double('z')

% z key:  Zoom in.

  zoom(2);

case double('Z')

% Z key:  Zoom out.

  zoom(0.5);

end

set(hfig,'UserData',handles);
if (moved)
  update_display(hfig);
end

return;
end  % KeyPressFcn

function change_type(hObject, eventdata, handles)
% hObject    handle to channel_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns channel_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from channel_list

hfig = gcbf;
contents = get(hObject,'String');
handles = get(hfig,'UserData');
handles.type = deblank(contents(get(hObject,'Value'),:));
set(hfig,'UserData',handles);

switch handles.type
  case 'markers'
    set(handles.cortex.pts(1), 'Marker','none','Color',[1 1 0]);
    %set(handles.cortex.pts(2), 'XData',NaN,'YData',NaN);
  otherwise
    set(handles.cortex.pts(1), 'Marker','o','Color',[1 0 0]);
end

update_display(hfig);

end

function update_display(hFig)

handles = get(hFig, 'UserData');
mymovie = handles.mymovie;
trackings = handles.trackings;

name = num2str(handles.current);

if (isempty(mymovie.(handles.type)))
  if (~isempty(mymovie.dic))
    handles.type = 'dic';
  elseif (~isempty(mymovie.markers))
    handles.type = 'markers';
  elseif (~isempty(mymovie.data))
    handles.type = 'data';
  end
  %ismember(get(handles.uitype, 'String'), handles.type)
  set(handles.uitype, 'Value', find(ismember(cellstr(get(handles.uitype, 'String')), handles.type)));
  %keyboard
  %set(handles.uitype, 'Value', handles.type);
end

if (isfield(mymovie, 'metadata'))
  %name = [name ':' num2str(mymovie.metadata.z_pos(handles.current)) ':' num2str(mymovie.metadata.z_rel(handles.current))];
end

set(hFig, 'Name', ['ASSET Results (' name ')']);

nans = NaN(1,2);

ref_egg = nans;
ref_cor = nans;
%err_egg = NaN;
%err_cor = NaN;
%det_egg = nans;
%det_cor = nans;
%pts_egg = nans;
%pts_cor = nans;
%ruf_cor = nans;
%ruf_prev = nans;
%ruf_int = nans;
%prev_pts = nans;
%prev_links = nans;

switch handles.type
%  case 'dic'
%    img = load_data(mymovie.dic, handles.current);
%    img2 = img;
%
%    if (isfield(mymovie, 'dic') & ~isempty(mymovie.dic))
%      pts_egg = mymovie.dic.centers(:,handles.current).';
%      pts_cor = [cos(mymovie.dic.orientations(:,handles.current)) sin(-mymovie.dic.orientations(:,handles.current))] * mymovie.dic.axes_length(1,handles.current) + pts_egg;
%      
%      if (isfield(mymovie.dic, 'ruffles'))
%        ruf_cor = mymovie.dic.ruffles(handles.current).carth;
%        ruf_prev = mymovie.dic.ruffles(handles.previous).carth;
%        ruf_indx = (all(mymovie.dic.ruffles(handles.current).bounds == 0, 2));
%        if (any(ruf_indx))
%          ruf_int = ruf_cor(ruf_indx,:);
%          ruf_cor = ruf_cor(~ruf_indx,:);
%        end
%      end
%    end
%
%    egg_type = 'dic';
%    cor_type = 'dic';

  case 'markers'

    fields = {'dic', 'dic', 'eggshell', 'cortex', 'eggshell', 'ruffles', 'ruffles'; ...
              'img', 'img', 'carth',    'carth',  'axes',     'carth',   'paths'};
    if (isfield(mymovie, 'eggshell') & ~isempty(mymovie.eggshell))
      fields{1,1} = 'eggshell';
    end
    if (isfield(mymovie, 'cortex') & ~isempty(mymovie.cortex))
      fields{1,2} = 'cortex';
    end

%    if (isfield(mymovie, 'eggshell') & ~isempty(mymovie.eggshell))
%      egg_type = 'markers';
%      img = load_data(mymovie.eggshell, handles.current);
%    else
%      egg_type = 'dic';
%      img = load_data(mymovie.dic, handles.current);
%    end
%    if (isfield(mymovie, 'cortex') & ~isempty(mymovie.cortex))
%      cor_type = 'markers';
%      img2 = load_data(mymovie.cortex, handles.current);
%    else
%      cor_type = 'dic';
%      img2 = load_data(mymovie.dic, handles.current);
%    end
%
%    if (isfield(mymovie, 'markers') & ~isempty(mymovie.markers))
%      pts_egg = mymovie.markers.centers(:,handles.current).';
%      pts_cor = [cos(mymovie.markers.orientations(:,handles.current)) sin(-mymovie.markers.orientations(:,handles.current))] * mymovie.markers.axes_length(1,handles.current) + pts_egg;
%
%      if (isfield(mymovie.markers, 'ruffles'))
%        ruf_cor = mymovie.markers.ruffles(handles.current).carth;
%        ruf_prev = mymovie.markers.ruffles(handles.previous).carth;
%        ruf_indx = all(mymovie.markers.ruffles(handles.current).bounds == 0, 2);
%        if (any(ruf_indx))
%          ruf_int = ruf_cor(ruf_indx,:);
%          ruf_cor = ruf_cor(~ruf_indx,:);
%        end
%      end
%    end
%  case 'data'
%    img = load_data(mymovie.data, handles.current);
%    img2 = img;
%
%    if (isfield(mymovie.data, 'centrosomes'))
%      pts_egg = mymovie.data.centrosomes(handles.current).carth(1,:);
%      pts_cor = mymovie.data.centrosomes(handles.current).carth(2,:);
%    end
%
%    egg_type = 'dic';
%    cor_type = 'dic';
  case 'data'
    fields = {'data', 'data', 'eggshell', 'cortex', 'eggshell', 'ruffles'; ...
              'img',  'img',  'carth',    'carth',  'axes',     'carth'};

    if (isfield(mymovie.data, 'spots') & ~empty_struct(mymovie.data.spots, 'carth'))
      fields(:,end+1) = {'spots'; 'carth'};
      if (~empty_struct(mymovie.data.spots, 'cluster'))
        fields(:,end+1) = {'spots'; 'cluster'};
      end
    end

  otherwise
    fields = {handles.type, handles.type, 'eggshell', 'cortex', 'eggshell', 'ruffles'; ...
              'img',        'img',        'carth',    'carth',  'axes',     'carth'};

end

if (~isempty(trackings) && isfield(trackings, handles.type))
  if (isfield(trackings.(handles.type), 'average'))
    ref_egg = trackings.(handles.type).average{1,handles.current};
    ref_cor = trackings.(handles.type).average{2,handles.current};
  elseif (isfield(trackings.(handles.type), 'eggshell') && isfield(trackings.(handles.type), 'cortex'))
    ref_egg = trackings.(handles.type).eggshell(handles.current).carth;
    ref_cor = trackings.(handles.type).cortex(handles.current).carth;
  end

  %tmp = mymovie.(handles.type).errors(1, handles.current, :, :);
  %err_egg = mymean(tmp(:));
  %tmp = mymovie.(handles.type).errors(2, handles.current, :, :);
  %err_cor = mymean(tmp(:));
end

values = extract_fields(mymovie, handles.current, handles.type, fields);

%det_egg = mymovie.(egg_type).eggshell(handles.current).carth;
%det_cor = mymovie.(cor_type).cortex(handles.current).carth;

%set(handles.eggshell.errors,'Data',{err_egg});
%set(handles.cortex.errors,'Data',{err_cor});

if (ishandle(handles.eggshell.img))
  set(handles.eggshell.img, 'CData', values{1});
  set(handles.cortex.img, 'CData', values{2});

  set(handles.eggshell.lines(1),'XData',ref_egg(:,1),'YData',ref_egg(:,2));
  set(handles.eggshell.lines(2),'XData',values{3}(:,1),'YData',values{3}(:,2));

  set(handles.cortex.lines(1),'XData',ref_cor(:,1),'YData',ref_cor(:,2));
  set(handles.cortex.lines(2),'XData',values{4}(:,1),'YData',values{4}(:,2));

  set(handles.eggshell.pts, 'XData', values{5}(:,1), 'YData', values{5}(:,2));

  set(handles.cortex.ruffles, 'XData', values{6}(:,1), 'YData', values{6}(:,2));

  if (numel(handles.cortex.pts) > 1 & handles.current ~= handles.previous)
    prev_pts = get(handles.cortex.pts(1), {'XData', 'YData'});

    if (~any(isnan(values{8})))
      goods = (values{8}(:,3) == handles.current - 1);
      if (handles.current > handles.previous)
        x = [values{7}(values{8}(goods,1),1), prev_pts{1}(values{8}(goods,2)).'];
        y = [values{7}(values{8}(goods,1),2), prev_pts{2}(values{8}(goods,2)).'];
      else
        prev_links = get(handles.cortex.pts(3), {'XData', 'YData'});
        if (numel(prev_links{1}) > 0)
          x = [values{7}(prev_links{1},1), prev_pts{1}(prev_links{2}).'];
          y = [values{7}(prev_links{1},2), prev_pts{2}(prev_links{2}).'];
        else
          x = [];
          y = [];
        end
      end
      
      x(:,end+1) = NaN;
      y(:,end+1) = NaN;

      x = x.';
      y = y.';

      set(handles.cortex.pts(2), 'XData', x(:), 'YData', y(:));

      goods = (values{8}(:,3) == handles.current - 1);
      set(handles.cortex.pts(3), 'XData', values{8}(goods,1), 'YData', values{8}(goods, 2));
    end
  end

  switch handles.type
    case 'markers'
      set(handles.cortex.pts(1), 'XData', values{7}(:,1), 'YData', values{7}(:,2));
    otherwise
      if (length(values) >= 7)
        set(handles.cortex.pts(1), 'XData', values{7}(:,1), 'YData', values{7}(:,2));
      end
%      set(handles.cortex.pts(2), 'XData', values{7}(:,1), 'YData', values{7}(:,2));
  end
  %set(handles.cortex.ruffles(2), 'XData', ruf_cor(:,1), 'YData', ruf_cor(:,2));
  %set(handles.cortex.ruffles(3), 'XData', ruf_int(:,1), 'YData', ruf_int(:,2));
else
  handles.eggshell.img = image(values{1},'Parent',handles.eggshell.axes,'CDataMapping','scaled');
  set(handles.eggshell.axes,'Visible','off');

  handles.cortex.img = image(values{2},'Parent',handles.cortex.axes,'CDataMapping','scaled');
  set(handles.cortex.axes,'Visible','off');

  handles.eggshell.lines(1,end+1) = line(ref_egg(:,1),ref_egg(:,2),'Color',[0 0 1],'Parent',handles.eggshell.axes);
  handles.eggshell.lines(1,end+1) = line(values{3}(:,1),values{3}(:,2),'Color',[0 1 0],'Parent',handles.eggshell.axes);

  handles.cortex.lines(1,end+1) = line(ref_cor(:,1),ref_cor(:,2),'Color',[1 0 0],'Parent',handles.cortex.axes);
  handles.cortex.lines(1,end+1) = line(values{4}(:,1),values{4}(:,2),'Color',[1 0.5 0],'Parent',handles.cortex.axes);

  handles.eggshell.pts = line(values{5}(:,1),values{5}(:,2),'Marker','o','Color',[1 1 0],'Parent',handles.eggshell.axes, 'LineStyle', 'none');

  handles.cortex.ruffles(1) = line(values{6}(:,1),values{6}(:,2),'Marker','*','Color',[0 0.5 0],'LineStyle','none','Parent',handles.cortex.axes);

  switch handles.type
    case 'markers'
      handles.cortex.pts(1) = line(values{7}(:,1),values{7}(:,2),'Color',[1 1 0],'Parent',handles.cortex.axes);
    case 'data'
      handles.cortex.pts(1) = line(values{7}(:,1),values{7}(:,2),'LineStyle','none','Marker','o','Color',[1 0 0],'Parent',handles.cortex.axes);
      if (numel(values) > 7)
        handles.cortex.pts(2) = line(NaN, NaN, 'Color', [0 0 1],'Parent',handles.cortex.axes);
        if (~any(isnan(values{8})))
          goods = ((values{8}(:,3) == handles.current - 1) & ~any(isnan(values{8}),2));
          handles.cortex.pts(3) = line(values{8}(goods,1), values{8}(goods,2), 'Visible', 'off','Parent',handles.cortex.axes);
        else
          handles.cortex.pts(3) = line(NaN, NaN, 'Visible', 'off','Parent',handles.cortex.axes);
        end
      end
    otherwise
      if (length(values) >= 7)
        handles.cortex.pts(1) = line(values{7}(:,1),values{7}(:,2),'Marker','o','Color',[1 0 0],'Parent',handles.cortex.axes);
      else
        handles.cortex.pts(1) = line(NaN, NaN,'Parent',handles.cortex.axes);
      end
%      handles.cortex.pts(2) = line(values{7}(2,1),values{7}(2,2),'Marker','o','Color',[0 0 1],'Parent',handles.cortex.axes);
  end
  %handles.cortex.ruffles(2) = line(ruf_cor(:,1),ruf_cor(:,2),'Marker','*','Color',[0 1 0],'LineStyle','none','Parent',handles.cortex.axes);
  %handles.cortex.ruffles(3) = line(ruf_int(:,1),ruf_int(:,2),'Marker','*','Color',[0 0 1],'LineStyle','none','Parent',handles.cortex.axes);
end

set(hFig,'UserData',handles);
end
