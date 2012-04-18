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
'Interruptible', 'off', ...
'Position',[10 10 200 60],...
'Visible','off');

hAxes = axes(...
'Parent',hFig,...
'Units','characters',...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'Interruptible', 'off', ...
'FontSize',10,...
'Tag','domain');
%'Position',[3.5 9.5 95 95],...

hFrameFig = figure(...
'Colormap',mygray, ....
'NumberTitle','off',...
'HandleVisibility','callback',...
'Tag','frame_fig',...
'UserData',[],...
'Visible','off');
%'MenuBar','none',...

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
  delete(hFrameFig);
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

case double('c')
% Clear the drawing

  setPosition(handles.polygon, NaN(1,2));
  %delete(handles.polygon);
  drawnow;
  [junk, x, y] = roipoly();

  setPosition(handles.polygon, [x,y]);
  uiwait(hfig);
  %handles.polygon = impoly(handles.axes, [x y], 'Closed', true);

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

case 32 % Space, confirm and next
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

  %disp(key)

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

  set(handles.display, 'Name', 'loading...');
  drawnow;

  name = handles.domains{handles.current};
  indx = strfind(name, 'DP');
  name = [name(1:indx(1)-2) '_.mat'];

  tmp = load(name);
  %img1 = imnorm(double(load_data(tmp.mymovie.dic, frame)));
  img2 = imnorm(double(load_data(tmp.mymovie.data, frame)));

  opts = get_struct('ASSET');
  [junk, pos] = gather_quantification(tmp.mymovie, opts);

  %pts = insert_ruffles(tmp.mymovie.markers.cortex(frame).carth, tmp.mymovie.markers.ruffles(frame).paths);
  pts = tmp.mymovie.data.quantification(frame).carth;
  [lin_pts, tot_length] = carth2linear(pts);
  [ell_pts] = carth2elliptic(pts, tmp.mymovie.markers.centers(:, frame), tmp.mymovie.markers.axes_length(:, frame), tmp.mymovie.markers.orientations(1, frame));

  center = lin_pts(find(abs(ell_pts(:,1)) == min(abs(ell_pts(:,1))), 1));
  lin_pts = lin_pts - center;
  lin_pts(lin_pts < -tot_length) = lin_pts(lin_pts < -tot_length) + tot_length;

  xy = getPosition(handles.polygon);
  xi = xy(:,1);
  yi = round(xy(:,2));

  x_pos = round(xi(yi == frame));

  if (~isempty(x_pos))

    tot_length = tot_length * opts.pixel_size;
    lin_pos = pos(x_pos);
    lin_pts = lin_pts * tot_length;

    dist = abs(cat(3, bsxfun(@minus, lin_pts, lin_pos), bsxfun(@minus, lin_pts - tot_length, lin_pos)));
    dist = min(dist, [], 3);
    [junk, indx] = min(dist, [], 1);

    pts = pts(indx, :);
    if (isempty(pts))
      pts = NaN(1, 2);
    end
  else
    pts = NaN(1,2);
  end

  childs = get(handles.display, 'Children');
  if (isempty(childs))
    %haxes = subplot(1,2,1, 'Parent', handles.display);
    %haxes(2) = subplot(1,2,2, 'Parent', handles.display);
    haxes = axes('Parent', handles.display);
    set(handles.display, 'Visible', 'on');
  else
    %haxes = childs(1:2);
    haxes = childs;
  end
  drawnow;
  %imshow(img1, 'Parent', haxes(1));
  imshow(img2, 'Parent', haxes);
  hold(haxes, 'on');
  scatter(haxes, pts(:, 1), pts(:, 2), 'r');
  hold(haxes, 'off');
  set(handles.display, 'Name', ['Frame ' num2str(frame)]);
  drawnow;

  %set(handles.display,'Visible', 'on');

  return;
end

function update_display(hFig)

handles = get(hFig, 'UserData');

if (handles.previous ~= handles.current)
  tmp = load(handles.domains{handles.current});
  handles.previous = handles.current;

  if (~isfield(tmp, 'domain'))
    if (isfield(tmp, 'img'))
      tmp.domain = tmp.img;
    elseif (isfield(tmp, 'mymovie'))
      tmp.domain = gather_quantification(tmp.mymovie, tmp.opts);
      tmp.path = NaN(0, 2);

      name = (handles.domains{handles.current});
      if (strncmp(name(end-2:end), 'mat', 3))
        name = name(1:end-4);
      end
      name = [name 'DP.mat'];

      if (exist(name) == 2)
        tmptmp = load(name);
        tmp.path = tmptmp.path;
      end

      handles.domains{handles.current} = name;

      tmp = rmfield(tmp, 'mymovie');
    end
  end

  if (~strncmp(class(tmp.domain), 'double', 6))
  img = double(tmp.domain);

  img(tmp.domain == intmax(class(tmp.domain))) = NaN;
  else
    img = tmp.domain;
  end
  clim = prctile(img(:), [0.1, 99.9]);
  img = imnorm(img, clim(1), clim(2));

  if (ishandle(handles.img))

    tmpx = get(handles.axes, 'XLim');
    tmpy = get(handles.axes, 'YLim');

    newx = [0 size(img, 2)];
    newy = [0 size(img, 1)];

    scaling = max(diff(newx)/diff(tmpx), diff(newy)/diff(tmpy));

    set(handles.img,'CData', img);
    zoom(1/scaling)

    set(handles.axes,'YLim',newy, 'XLim',newx);
    handles.zoom = 1;
  else

    %handles.img = imagesc(img,'Parent', handles.axes,'CDataMapping', 'scaled');
    handles.img = imagesc(img,'Parent', handles.axes, 'ButtonDownFcn', @mouseClick);
    %aspect(handles.axes, [aspect_ratio]);
    set(handles.axes,'Visible', 'off');
    %colormap(handles.axes, jet);
    j = jet(128);
    set(hFig, 'Colormap', j);
  end
  set(hFig,'Visible', 'on', 'Name', handles.domains{handles.current});

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

  drawnow
%  drawnow
set(hFig,'UserData',handles);
end

end

function store_path(hFig)

  %keyboard

  handles = get(hFig, 'UserData');

  domain = get(handles.img, 'CData');
  [h,w] = size(domain);

  xy = getPosition(handles.polygon);

  if (all(isnan(xy(:))))
    btn = questdlg('No segmentation detected, what do you want to do ? Note that you need to double-click on the polygon to save it...', 'Manual segmentation ?', 'Ignore', 'Cancel', 'Cancel');

    switch btn
      case 'Ignore'
        path = NaN(1,2);
        save(handles.domains{handles.current}, 'domain', 'path');
      case 'Cancel'
        return;
    end
  else

  [tmp_xy, i, j] = unique(xy, 'rows');

  xy = NaN(size(xy));
  xy(i, :) = tmp_xy;
  xy = xy(~any(isnan(xy), 2), :);

  xi = xy(:,1);
  yi = floor(xy(:,2));

  yi(yi < 1) = 1;
  yi(yi > h) = h;

    start = find(yi == min(yi), 1);
    if (start ~= 1)
      xi = [xi(start:end); xi(1:start-1)];
      yi = [yi(start:end); yi(1:start-1)];
    end
    if (yi(1) == yi(2) & xi(1) == xi(2))
      xi = xi(2:end);
      yi = yi(2:end);
    end
    if (yi(1) >= yi(2))
      xi = xi([1 end:-1:2]);
      yi = yi([1 end:-1:2]);
    end
    start = yi(1);

    dx = abs(diff(xi));
    ends = find(yi == max(yi), 1);
    if (dx(ends) < dx(ends-1))
      ends = ends - 1;
    end

    %ends = find(dy ~= 1, 1, 'first');
  
    half1 = [xi(1:ends) yi(1:ends)];
    half2 = [xi([ends+1:end 1]) yi([ends+1:end 1])];
    half2 = half2([end:-1:1], :);
    if (half2(1, 2) ~= half1(1, 2))
      half2 = [half1(1,:); half2];
    end
    half1(end, 2) = h;
    half2(end, 2) = h;

    [junk, indx] = unique(half1(:, 2));
    half1 = half1(indx, :);
    [junk, indx] = unique(half2(:, 2));
    half2 = half2(indx, :);

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
    
    %figure;scatter(path(:, 1), path(:, 2))

    save(handles.domains{handles.current}, 'domain', 'path');
  end


  return;
end
