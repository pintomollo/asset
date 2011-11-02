function manually_segment_domains(fname)

  if (ischar(fname))
    if (~isempty(findstr(fname, '*')))
      tmp = dir(fname); 
      mymovies = cell(1, length(tmp));

      for i=1:length(tmp)
        mymovies{i} = tmp(i).name;
      end
    else
      mymovies = {fname}; 
    end
  elseif (~iscell(fname))
    mymovies = {};
    mymovies{1} = fname;
  end

mygray = [0:255]' / 255;
mygray = [mygray mygray mygray];

hFig = figure(...
'Units','characters',...
'Color',[0.701960784313725 0.701960784313725 0.701960784313725],...
'Colormap',mygray, ....
'KeyPressFcn',@KeyPressFcn, ...
'CloseRequestFcn',@FexitFcn, ...
'MenuBar','none',...
'NumberTitle','off',...
'HandleVisibility','callback',...
'Tag','domains_fig',...
'Resize','off',...
'UserData',[],...
'Visible','off');
%'Position',[104 6 212 54],...

hAxes = axes(...
'Parent',hFig,...
'Units','characters',...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'FontSize',10,...
'Tag','domain');
%'Position',[3.5 9.5 95 95],...

hFrameFig = figure(...
'Colormap',mygray, ....
'MenuBar','none',...
'NumberTitle','off',...
'HandleVisibility','callback',...
'Tag','frame_fig',...
'UserData',[],...
'Visible','off');

  handles = struct('display', hFrameFig,  ...
                   'axes', hAxes, ...
                   'domains', {mymovies}, ...
                   'polygon', -1, ...
                   'img', -1, ...
                   'previous', 0, ...
                   'zoom', 1, ...
                   'current', 1);

  set(hFig, 'UserData',  handles);

  update_display(hFig);

  uiwait(hFig);
  delete(hFig);
  drawnow;

  return;
end
%  nmovies = length(mymovies);

%  for i=1:nmovies
%    mymovie = mymovies{i};
%
%    if (ischar(mymovie))
%      domain = imread(mymovie);
%    elseif (isnumeric(mymovie))
%      domain = mymovie;
%      mymovie = 'unkown_name.png'
%    end
%    [h, w] = size(domain);

%    imagesc(domain);
%    [bw, xi, yi] = roipoly();
%
%    if (numel(xi) < 10)
%      continue;
%    end
%
%    start = find(yi == min(yi), 1);
%    if (start ~= 1)
%      xi = [xi(start:end); xi(1:start-1)];
%      yi = [yi(start:end); yi(1:start-1)];
%    end
%    start = ceil(yi(1));
%
%    dy = [0; sign(diff(yi))];
%    ends = find(dy == -dy(2), 1, 'first');
%    dist = hypot(diff(xi(ends-1:ends+1)), diff(yi(ends-1:ends+1)));
%  
%    if (dist(1) > dist(2))
%      ends = ends-1;
%    end
%
%    half1 = [xi(1:ends) yi(1:ends)];
%    half2 = [xi([ends+1:end 1]) yi([ends+1:end 1])];
%
%    if (dy(2) < 0)
%      half1 = half1([end:-1:1], :);
%    else
%      half2 = half2([end:-1:1], :);
%    end
%
%    half1(end, 2) = h;
%    half2(end, 2) = h;
%
%    indexes = [start:h].';
%    x1 = interp1q(half1(:, 2), half1(:, 1), indexes);
%
%    x2 = interp1q(half2(:, 2), half2(:, 1), indexes);
%
%    mean_pos = mean([x1,x2], 2);
%    width = abs((x1 - x2) / 2);
%
%    path = NaN(h, 2);
%    path(start:end,:) = [mean_pos, width];
%
%    %imagesc(domain)
%    %hold on;
%    %plot(path(:, 1), [1:h], 'k');
%    %plot(path(:, 1) + path(:,2), [1:h], 'k');
%    %plot(path(:, 1) - path(:,2), [1:h], 'k');
%    %hold off;
%
    %keyboard
%
%    save([mymovie(1:end-3) '-Manual-DP.mat'], 'domain', 'path');
%  end
%
%  return;
%end

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

hfig = ancestor(hobj, 'figure');
handles = get(hfig,'UserData');

switch key
  case 28

    if (~isempty(modifier))
      handles.previous = handles.current;

      handles.current = handles.current - 1;
      if (handles.current < 1)
        handles.current = length(handles.domains);
      end
    else

      xaxes = get(handles.axes,'XLim');
      xaxes = xaxes - (diff(xaxes) /4);

      xaxes = xaxes - (xaxes(1) * (xaxes(1) < 0));

      set(handles.axes,'XLim',xaxes);
    end

% Left Arrow keypress:  rotate camera clockwise (as viewed
% from positive z).

  case 29

% Right Arrow keypress:  rotate camera ccw.

  if (~isempty(modifier))
      keyboard

      handles.previous = handles.current;
      handles.current = handles.current + 1;
      if (handles.current > length(handles.domains))
        handles.current = 1;
      end
  else
      maxx = size(get(handles.img,'CData'),2);

      xaxes = get(handles.axes,'XLim');
      xaxes = xaxes + (diff(xaxes) /4);

      xaxes = (xaxes - (xaxes(2) - maxx) * (xaxes(2) > maxx)) ;

      set(handles.axes,'XLim',xaxes);
    end

case 30

% Up Arrow keypress:  increase elevation.
  
  if (isempty(modifier))

    yaxes = get(handles.axes,'YLim');
    yaxes = yaxes - (diff(yaxes) /4);
      
    yaxes = yaxes - (yaxes(1) * (yaxes(1) < 0));

    set(handles.axes,'YLim',yaxes);
  end

case 31

% Down Arrow keypress:  decrease elevation.
  if (isempty(modifier))
    maxy = size(get(handles.img,'CData'),1);

    yaxes = get(handles.axes,'YLim');
    yaxes = yaxes + (diff(yaxes) /4);

    yaxes = (yaxes - (yaxes(2) - maxy) * (yaxes(2) > maxy));

    set(handles.axes,'YLim',yaxes);
  end

case double('z')

% z key:  Zoom in.

  zoom(handles.axes, 2);
  handles.zoom = handles.zoom * 2;

case double('Z')

% Z key:  Zoom out.

  zoom(handles.axes,0.5);
  handles.zoom = handles.zoom * 0.5;
  if (handles.zoom < 1)
    handles.zoom = 1;
  end

case 13 % Enter, confirm and next
      store_path(hfig);

      handles.previous = handles.current;
      handles.current = handles.current + 1;
      if (handles.current > length(handles.domains))
        handles.current = 1;
      end

case 27 % ESC, exit

  FexitFcn(hfig, []);


case 49 % 1, zoom normal

  zoom(handles.axes, 1/handles.zoom);
  handles.zoom = 1;

otherwise

  disp(key)

end

set(hfig,'UserData',handles);
update_display(hfig);

return;
end  % KeyPressFcn

function FexitFcn(hfig, event_data)

  store_path(hfig);
  uiresume(hfig);

  return;
end

function mouseClick(hobj, event_data)

  %handles = get(hFig, 'UserData');
  %get(handles.axes, 'CurrentPoint')
  hfig = ancestor(hobj, 'figure');
  pos = get(get(hobj, 'Parent'), 'CurrentPoint');
  frame = floor(pos(1,2));
  if (frame < 1)
    frame = 1;
  end

  handles = get(hfig, 'UserData');
  name = handles.domains{handles.current};
  indx = strfind(name, '.');
  name = name(1:indx(1)-1);
  indx = strfind(name, '-');
  name = [name(1:indx(end)-1) '_.mat'];

  tmp = load(name);
  img1 = imnorm(double(load_data(tmp.mymovie.dic, frame)));
  img2 = imnorm(double(load_data(tmp.mymovie.data, frame)));

  childs = get(handles.display);
  if (isempty(childs))
    haxes = subplot(121, 'Parent', handles.display);
    haxes(2) = subplot(122, 'Parent', handles.display);
    %haxes = axes('Parent', handles.display);
    set(handles.display, 'Visible', 'on');
  else
    haxes = childs(1:2);
  end
  drawnow;
  imshow(img1, 'Parent', haxes(1));
  imshow(img2, 'Parent', haxes(2));

  return;
end

function update_display(hFig)

handles = get(hFig, 'UserData');

if (handles.previous ~= handles.current)
  tmp = load(handles.domains{handles.current});
  handles.previous = handles.current;

  if (~isfield(tmp, 'domain'))
    tmp.domain = tmp.img;
  end

  if (~strncmp(class(tmp.domain), 'double', 6))
  img = double(tmp.domain);

  img(tmp.domain == intmax(class(tmp.domain))) = NaN;
  else
    img = tmp.domain;
  end
  img = imnorm(img);

  if (ishandle(handles.img))
    set(handles.img,'CData', img);
  else

    %handles.img = imagesc(img,'Parent', handles.axes,'CDataMapping', 'scaled');
    handles.img = imagesc(img,'Parent', handles.axes, 'ButtonDownFcn', @mouseClick);
    %aspect(handles.axes, [aspect_ratio]);
    set(handles.axes,'Visible', 'off');
    %colormap(handles.axes, jet);
    set(hFig, 'Colormap', jet);
  end
  set(hFig,'Visible', 'on');

  pos = tmp.path(:,1) - tmp.path(:,2);
  pos = [pos, [1:length(pos)].'];
  rev_pos = tmp.path(:,1) + tmp.path(:,2);
  rev_pos = [rev_pos, [1:length(rev_pos)].'];

  pos = [pos; rev_pos(end:-1:1,:)];
  pos = pos(~(any(isnan(pos), 2)), :);
  if (handles.polygon ~= -1)
    setPosition(handles.polygon, pos);
  else
  handles.polygon = impoly(handles.axes, pos, 'Closed', true);
  end

  refresh
  drawnow
    colormap(jet);
%  drawnow
set(hFig,'UserData',handles);
end

end

function store_path(hFig)

  handles = get(hFig, 'UserData');

  domain = get(handles.img, 'CData');
  [h,w] = size(domain);

  xy = getPosition(handles.polygon);
  xi = xy(:,1);
  yi = xy(:,2);

     start = find(yi == min(yi), 1);
    if (start ~= 1)
      xi = [xi(start:end); xi(1:start-1)];
      yi = [yi(start:end); yi(1:start-1)];
    end
    start = round(yi(1));

    dy = [0; sign(diff(yi))];
    ends = find(dy == -dy(2), 1, 'first');
    dist = hypot(diff(xi(ends-1:ends+1)), diff(yi(ends-1:ends+1)));
  
    if (dist(1) > dist(2))
      ends = ends-1;
    end

    half1 = [xi(1:ends) yi(1:ends)];
    half2 = [xi([ends+1:end 1]) yi([ends+1:end 1])];

    if (dy(2) < 0)
      half1 = half1([end:-1:1], :);
    else
      half2 = half2([end:-1:1], :);
    end
    half1(end, 2) = h;
    half2(end, 2) = h;

    indexes = [start:h].';
    x1 = interp1q(half1(:, 2), half1(:, 1), indexes);

    x2 = interp1q(half2(:, 2), half2(:, 1), indexes);

    mean_pos = mean([x1,x2], 2);
    width = abs((x1 - x2) / 2);

    path = NaN(h, 2);
    path(start:end,:) = [mean_pos, width];

    %tmp = load(handles.domains{handles.current});
    %[tmp.path path]
    %keyboard
    save(handles.domains{handles.current}, 'domain', 'path');

  return;
end
