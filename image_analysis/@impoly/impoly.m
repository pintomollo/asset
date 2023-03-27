function h = impoly(varargin)
% IMPOLY draws a draggable and editable polygon. Freely inspired by the MATLAB
% class of the same name, to allow code intercompatibility. Functionalities
% should be similar.
%
%   P = IMPOLY(HAX, POS)    Draws an editable and draggable polygon on axes
%                           HAX using the Nx2 coordinates POS and returns an IMPOLY object P. 
%   IMPOLY(POS)             When no handle is provided, the current axes are used.
%   IMPOLY(HAX)             When the coordinates are not provided, the polygon is
%                           interactively drawn.
%   IMPOLY()                Interactive drawing on the current axes.
%   IMPOLY(..., 'Closed', BOOL) % Determines whether the polygon is a closed shape or not (FALSE is default).
%   IMPOLY(..., 'PositionConstraintFcn', FUNC) % Not currently implemented.
%   IMPOLY(..., 'Color', COL)   % Determines the color of the polygon ('k' is default).
%
% Interactions with the polygon:
%   Creation mode:
%       LEFT-CLICK          : Add a new vertex to the polygon
%       RIGHT-CLICK         : Toggles whether the poylgon is closed or not
%       MIDDLE-CLICK        : Exit creation mode
%       DOUBLE LEFT-CLICK   : Adds a last vertex and exit creation mode
%   Edition mode:
%       LEFT-CLICK          : If on a vertex, moves that vertex
%                           : If on an edge, drags the polygon
%       RIGHT-CLICK         : If on a vertex, deletes that vertex
%                           : If on an edge, creates a new vertex there
%
% Simon Blanchoud
% University of Fribourg
% 01.10.18

  % Check the input parameters
  [reg, is_closed, constraint, color] = parseparams(varargin, 'Closed', ...
        false, 'PositionConstraintFcn', [], 'Color', 'k');

  % Sets the default values for the axes and position
  if (isempty(reg))
    hax = gca();
    pos = [];
  elseif (numel(reg) == 1)
    if (ishghandle(reg{1}))
      % Always look for an axis somewhere in the ancestors list
      hax = ancestor(reg{1}, 'axes');
      pos = [];
    else
      hax = gca ();
      pos = reg{1};
    end
  else
    if (numel(reg) > 2)
      error('impoly:invalidInputs', 'Invalid number of arguments.');
    end

    hax = ancestor(reg{1}, 'axes');
    pos = reg{2};
  end

  % Require a valid handle to an axis
  if (isempty(hax))
    error('impoly:invalidInputs', 'Not a valid axes handle.')
  end

  % Prepare the interactive creation mode
  is_interactive = false;
  if (isempty(pos))
    pos = NaN(1,2);
    is_interactive = true;

  % We don't want duplicated end points, impoly handles that internally
  elseif (all(pos(1,:) == pos(end,:)))
    pos = pos(1:end-1, :);
  end

  % Create the polygon, with vertexes and an additional property to check for completion
  hlin = line(pos(:,1), pos(:,2), 'parent', hax, 'color', color, 'marker', '+');
  addproperty('is_done', hlin, 'boolean', false);

  % The parameter structure used to store the information on the polygon. The trick here
  % is that the object 'impoly' is basically a handle to the graphical object, in which
  % the actual data is store. Which means that the 'impoly' object is interactively updated too.
  mypoly = struct('hPolygon', hlin, ...
                  'is_closed', is_closed, ...
                  'color', color, ...
                  'position', pos, ...
                  'init', [], ...
                  'is_draggable', true);

  % Store that structure in the graphical object and convert it to impoly for returning
  set(hlin, 'userdata', mypoly);
  mypoly = class(mypoly, 'impoly');

  % Enter interactive creation mode
  if (is_interactive)

    % Empty all callbacks to prevent other functions to run in parallel
    [hfig, bkp] = lock_window(hlin, hax, mypoly);
    set(hfig, 'WindowButtonDownFcn', @create_polygon);
    set(hfig, 'WindowButtonMotionFcn', @draw_polygon);

    % Wait until the polygon is created
    waitfor(hlin, 'is_done');
    
    % Restore the callbacks
    lock_window(hfig, bkp);
  end

  % Set the interactive edition of the polygon
  set(hlin, 'buttondownfcn', @select_polygon);

  % Return the impoly object 
  if (nargout > 0)
    h = mypoly;
  end

  return;
end

function create_polygon(hsrc, evt)
% Callback used to create the polygon with each consecutive mouse click

  % Get the data of the polygon hidden in the graphic object
  data = get(hsrc, 'userdata');
  clickType = get(hsrc, 'SelectionType');

  % Adapt the behavior to the type of mouse click
  switch clickType

    % Left-click, add a point to the polygon
    case 'normal'

      % Get the current position of the mouse
      pos = get(data.axes, 'currentpoint');

      % And the vertices of the polygon
      vert = getPosition(data.impoly);

      % This represents an uninitialized polygon
      if (any(isnan(vert(:))))
        vert = pos(1,1:2);
      else
        vert = [vert; pos(1,1:2)];
      end

      % Store the new vertices
      setPosition(data.impoly, vert);

    % Right-click toggles the closed/open state of the polygon
    case 'alt'
      setClosed(data.impoly, ~getClosed(data.impoly));

      % Updates the polygon to the current status of edition
      draw_polygon(hsrc);

    % Middle-click or double-click finish the polygon
    case  {'extend', 'open'}

      % Notify the end and refresh the drawn polygon
      set(data.hPolygon, 'is_done', true);
      setPosition(data.impoly, getPosition(data.impoly));
  end

  return;
end

function draw_polygon(hsrc, evnt)  %#ok
% Callback function used to draw the current polygon, meaning the one the
% user is drawing. So that's the stored one plus the next point to be added.
            
  % Get the data stored for the drawing in the figure
  data = get(hsrc, 'userdata');
  % Get the internal data of the polygon
  poly = get(data.hPolygon, 'userdata');
  % And the current position of the mouse
  pos = get(data.axes, 'currentpoint');

  if poly.is_closed
      xpos = [poly.position(:, 1); pos(1,1); poly.position(1, 1)];
      ypos = [poly.position(:, 2); pos(1,2); poly.position(1, 2)];
  else
      xpos = [poly.position(:, 1); pos(1,1)];
      ypos = [poly.position(:, 2); pos(1,2)];
  end

  set(data.hPolygon, 'xdata', xpos, 'ydata', ypos);
  drawnow();

  return;
end

function select_polygon(hsrc, evt)

  data = get(hsrc, 'userdata');
  haxes = ancestor(hsrc, 'axes');
  pos = get(haxes, 'currentpoint');
  prev = data.position;

  u = get(haxes, 'units');
  set(haxes, 'units', 'pixels');
  dp = get(haxes, 'position');
  set(haxes, 'units', u);
  dm = get(data.hPolygon, 'MarkerSize')*4/3;
  dv = axis(haxes);

  prange = dp(3:4) - dp(1:2);
  vrange = dv([2 4]) - dv([1 3]);
  scaling = max(prange./vrange);

  x = pos(1,1);
  y = pos(1,2);
  xpos = prev(:,1);
  ypos = prev(:,2);

  [vdist, vind] = min((xpos - pos(1,1)) .* (xpos - pos(1,1)) + ...
  (ypos - pos(1,2)) .* (ypos - pos(1,2)));

  if (vdist*scaling*scaling < dm*dm)

    if (evt == 1)
    pause(0.01);
    drawnow();
    pause(0.01);

      [hfig, bkp] = lock_window(hsrc, haxes);

      data.init = vind;
      set(hsrc, 'userdata', data);

      set(hfig, 'WindowButtonMotionFcn', @edit_polygon, ...
          'WindowButtonDownFcn', @waitforclick);

      %[x, y, btn] = ginput(1)
      waitfor(hsrc, 'is_done');
      lock_window(hfig, bkp);

      data = get(hsrc, 'userdata');
      %x = data.init(1);
      %y = data.init(2);

      if (~isempty(data.init))
      %if (btn == 1)
          xpos(vind) = data.init(1);
          ypos(vind) = data.init(2);

      end
        update_polygon(data, xpos, ypos);

      data.init = [];
    elseif (evt == 3)

      xpos(vind) = [];
      ypos(vind) = [];

      if (isempty(xpos))
        delete(hsrc);
      else
        update_polygon(data, xpos, ypos);
      end
    end

  else

    if (evt == 1)
    pause(0.01);
    drawnow();
    pause(0.01);
      [hfig, bkp] = lock_window(hsrc, haxes);

      data.init = pos(1,1:2);
      set(hsrc, 'userdata', data);

      set(hfig, 'WindowButtonMotionFcn', @edit_polygon, ...
          'WindowButtonDownFcn', @waitforclick);


      waitfor(hsrc, 'is_done');
      %[x, y, btn] = ginput(1)
      lock_window(hfig, bkp);
      data = get(hsrc, 'userdata');

      if (~isempty(data.init))
        xpos = xpos + data.init(1) - pos(1,1);
        ypos = ypos + data.init(2) - pos(1,2);
      end

      update_polygon(data, xpos, ypos);
    elseif (evt == 3)

      if (data.is_closed)
        xpos = xpos([1:end 1]);
        ypos = ypos([1:end 1]);
      end

      xvect = xpos(2:end) - xpos(1:end-1);
      yvect = ypos(2:end) - ypos(1:end-1);

      [dist, indx] = min(abs(yvect*x - xvect*y + xpos(2:end).*ypos(1:end-1) - ...
      ypos(2:end).*xpos(1:end-1)) ./ ...
      (yvect.*yvect + xvect.*xvect));
      xpos = [xpos(1:indx); x; xpos(indx+1:end)];
      ypos = [ypos(1:indx); y; ypos(indx+1:end)];

      if (data.is_closed)
        xpos = xpos(1:end-1);
        ypos = ypos(1:end-1);
      end

      update_polygon(data, xpos, ypos);
    end
  end                  

  return;
end

function edit_polygon(hsrc, evnt)  %#ok
            
    data = get(hsrc, 'userdata');
    poly = get(data.hPolygon, 'userdata');
    pos = get(data.axes, 'currentpoint');
    
    xpos = poly.position(:, 1);
    ypos = poly.position(:, 2);
    if (numel(poly.init)==1)
        xpos(poly.init) = pos(1,1);
        ypos(poly.init) = pos(1,2);
    else
        xpos = xpos + (pos(1,1) - poly.init(1));
        ypos = ypos + (pos(1,2) - poly.init(2));
    end

    if poly.is_closed
        xpos = xpos([1:end 1]);
        ypos = ypos([1:end 1]);
    end

    set(data.hPolygon, 'xdata', xpos, 'ydata', ypos);
    refresh(hsrc);
end


function waitforclick(hsrc, evt)

data = get(hsrc, 'userdata');
poly = get(data.hPolygon, 'userdata');
pos = get(data.axes, 'currentpoint');

if (evt == 1)
  poly.init = pos(1, 1:2);
else
  poly.init = [];
end
set(data.hPolygon, 'userdata', poly, 'is_done', true);

  return;
end

function update_polygon(data, xpos, ypos)

  data.position = [xpos(:), ypos(:)];

  if (data.is_closed)
      xpos = xpos([1:end 1]);
      ypos = ypos([1:end 1]);
  end
  set(data.hPolygon, 'xdata', xpos, 'ydata', ypos);

  set(data.hPolygon, 'userdata', data);

  return;
end

function [hfig, bkp] = lock_window(hpoly, haxes, mypoly)

  % Unlock the window
  if (nargout == 0)
    hfig = hpoly;
    bkp = haxes;

    set(hfig, 'userdata', bkp.bkp, ...
        'ButtonDownFcn', bkp.dwnfnc, ...
        'KeyPressFcn', bkp.kpressfnc, ...
        'KeyReleaseFcn', bkp.kreleafnc, ...
        'SizeChangedFcn', bkp.szchgfnc, ...
        'WindowButtonDownFcn', bkp.wbtnfnc, ...
        'WindowButtonMotionFcn', bkp.wmovfnc, ...
        'WindowButtonUpFcn', bkp.wupfnc, ...
        'WindowKeyPressFcn', bkp.wkpressfnc, ...
        'WindowKeyReleaseFcn', bkp.wkreleafnc, ...
        'WindowScrollWheelFcn', bkp.wscrllwfnc); 

    set(bkp.hon, 'hittest', 'on');

  else
    if (nargin < 3)
      mypoly = [];
    end

    % Gets all the interactive objects to we can turn them off temporarily
    hon = findobj('hittest', 'on');
    %hon = hon(hon ~= hpoly && hon ~= haxes);
    set(hon, 'hittest', 'off');

    % We also need the figure to replace its various interaction handles
    hfig = ancestor(haxes, 'figure');

    % Save everything we will be replacing
    bkp = struct('dwnfnc', get(hfig, 'ButtonDownFcn'), ...
        'kpressfnc',  get(hfig, 'KeyPressFcn'), ...
        'kreleafnc',  get(hfig, 'KeyReleaseFcn'), ...
        'szchgfnc',   get(hfig, 'SizeChangedFcn'), ...
        'wupfnc',     get(hfig, 'WindowButtonUpFcn'), ...
        'wmovfnc',    get(hfig, 'WindowButtonMotionFcn'), ...
        'wbtnfnc',    get(hfig, 'WindowButtonDownFcn'), ...
        'wkpressfnc', get(hfig, 'WindowKeyPressFcn'), ...
        'wkreleafnc', get(hfig, 'WindowKeyReleaseFcn'), ...
        'wscrllwfnc', get(hfig, 'WindowScrollWheelFcn'), ...
        'axes', haxes, ...
        'bkp', get(hfig, 'userdata'), ...
        'hon', hon, ...
        'hPolygon', hpoly, ...
        'impoly', mypoly);

    set(hfig, 'userdata', bkp, ...
        'ButtonDownFcn', [], ...
        'KeyPressFcn', [], ...
        'SizeChangedFcn', [], ...
        'WindowButtonDownFcn', [], ...
        'WindowButtonMotionFcn', [], ...
        'WindowKeyPressFcn', [], ...
        'WindowScrollWheelFcn', []);

    %    'KeyReleaseFcn', [], ...
    %    'WindowKeyReleaseFcn', [], ...
    %    'WindowButtonUpFcn', [], ...
    %refresh(hfig);

    set(hpoly, 'is_done', false);
  end

  return;
end

