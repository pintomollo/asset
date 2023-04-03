function dragzoom2D(hAx)
%DRAGZOOM2D Drag and zoom tool (simplified from dragzoom.m by Evgeny Pr)
%
% Description:
%   DRAGZOOM2D allows you to interactively manage the axes in a figure
%   through draging and zooming using the mouse and keyboard shortcuts.
%
% Using:
%   dragzoom()                      : manages the axis from gca()
%   dragzoom(hAx)                   : manages the axis from the handle hAx
%
% Interactions:
%   Mouse actions:
%       single-click and holding LB : Activation Drag mode
%       single-click and holding RB : Activation Rubber Band for region zooming
%       single-click MB             : Activation 'Extend' Zoom mode
%       double-click LB, RB, MB     : Reset to Original View
% 
%   Hotkeys:
%       '+' or '>' or ','           : Zoom plus
%       '-' or '<' or '.'           : Zoom minus
%       '0' or '='                  : Set default axes (reset to original view)
%       'uparrow'                   : Up drag
%       'downarrow'                 : Down drag
%       'leftarrow'                 : Left drag
%       'rightarrow'                : Right drag
%       'x'                         : Toggles zoom and drag only for X axis
%       'y'                         : Toggles zoom and drag only for Y axis
%
% Adapted from dragzoom.m, version 0.9.7 by Evgeny Pr aka iroln <esp.home@gmail.com>
% Simon Blanchoud
% University of Fribourg
% 28.09.2018

% The initial objective was to make dragzoom.m compatible with Octave by solving the 
% nested callback problems. Given the tremendous amount of code written by Evgeny, I
% decided to simplify it as much as possible and stick to 2D. In addition, I tried
% to compact the original code and linearize it as much as possible. Original calls
% to linearized functions are indicated by "%---". To circumvent the limitations of
% removing nested functions, a parameter structure has been created and is stored in
% the parent figure.

    % Get the handles to the axis and the parent figure
    if (nargin == 0)
        hAx = gca();
        hFig = ancestor(hAx, 'figure');
    elseif (nargin == 1 && ishandle(hAx) && strncmp(get(hAx, 'Type'), 'axes', 4) && length(hAx) == 1)
        hFig = ancestor(hAx, 'figure');
    else
        error('dragzoom2D:invalidInputs', ...
            'Input must be a single axes handle.');
    end

    %---params = Setup(hAx);
    % Setup options
    % Used to test whether the axis has an image in it
    h = findobj(hAx, 'Type', 'Image');
        
    % Default values for the zoom grid
    mZoomMinPow = 0;
    mZoomMaxPow = 5;
    mZoomNum = 51;
    mZoomIndex = 11; % index of 100%
    [mDefaultZoomGrid, mDefaultZoomSteps] = ...
        ZoomLogGrid(mZoomMinPow, mZoomMaxPow, mZoomNum);
    
    % Parameters structure with handles and default values
    params = struct('hAx', hAx, ...
        'mStartX', [], ...
        'mStartY', [], ...
        'mBindX', [], ...
        'mBindY', [], ...
        'mDragShiftStep', 3, ...
        'mDragSaveShiftStep', 3, ...
        'mDragShiftStepInc', 1, ...
        'mZoomMinPow', mZoomMinPow, ...
        'mZoomMaxPow', mZoomMaxPow, ...
        'mZoomNum', mZoomNum, ...
        'mZoomExtendNum', 301, ...
        'mZoomKeysNum', 181, ...
        'mDefaultZoomGrid', mDefaultZoomGrid, ...
        'mDefaultZoomSteps', mDefaultZoomSteps, ...
        'mZoomGrid', mDefaultZoomGrid, ...
        'mZoomSteps', mDefaultZoomSteps, ...
        'mZoomIndexX', mZoomIndex, ...
        'mZoomIndexY', mZoomIndex, ...
        'mDefaultXLim', get(hAx, 'XLim'), ...
        'mDefaultYLim', get(hAx, 'YLim'), ...
        'mRubberBand', [], ...
        'mRbEdgeColor', 'k', ...
        'mRbFaceColor', 'none', ...
        'mRbFaceAlpha', 1, ...
        'fIsDragAllowed', false, ...
        'fIsZoomExtendAllowed', false, ...
        'fIsRubberBandOn', false, ...
        'fIsEnableDragX', true, ...
        'fIsEnableDragY', true, ...
        'fIsEnableZoomX', true, ...
        'fIsEnableZoomY', true, ...
        'fIsImage', ~isempty(h));

    % Defines the handles to the callback functions and the parameter structure
    set(hFig, 'CurrentAxes', hAx, ...
        'userdata', params, ...
        'WindowButtonDownFcn',      {@WindowButtonDownCallback2D}, ...
        'WindowButtonUpFcn',        {@WindowButtonUpCallback2D}, ...
        'WindowButtonMotionFcn',    {@WindowButtonMotionCallback2D}, ...
        'KeyPressFcn',              {@WindowKeyPressCallback2D}, ...
        'KeyReleaseFcn',            {@WindowKeyReleaseCallback2D}, ...
        'WindowKeyPressFcn',        {@WindowKeyPressCallback2D}, ...
        'WindowKeyReleaseFcn',      {@WindowKeyReleaseCallback2D});

    return;
end
%--------------------------------------------------------------------------

%==========================================================================
function WindowButtonDownCallback2D(src, evnt)    %#ok
    %WindowButtonDownCallback2D called when the mouse clicks. Typically,
    % mouse position will be stored to be compared during the movement to 
    % update the axis according to the amplitude of the movement.
    
    % Retrieve the parameter structure and the type of event
    params = get(src, 'userdata');
    clickType = get(src, 'SelectionType');
    
    % Defines the type of mouse click
    switch clickType

        % Left click: dragging
        %case 'normal'
        %    %----DragMouseBegin();
        %    if (~params.fIsDragAllowed)
        %        [cx, cy] = GetCursorCoordOnWindow(src);
        %        
        %        params.mStartX = cx;
        %        params.mStartY = cy;
        %        
        %        params.fIsDragAllowed = true;
        %    end

        % Double click: reset
        case 'open'
            %---ResetAxesToOrigView();
            SetAxesLimits(params.hAx, params.mDefaultXLim, params.mDefaultYLim);
            params.mZoomIndexX = find(params.mZoomGrid == 100);
            params.mZoomIndexY = params.mZoomIndexX;

        % Right click: rubber band zoom
        case 'alt'
            %---RubberBandBegin();
            if (~params.fIsRubberBandOn)
                [acx, acy] = GetCursorCoordOnAxes(params.hAx);
               
                %---RubberBandSetup();
                h1 = patch([acx acx acx acx], [acy acy acy acy], 'k',
                    'Parent', params.hAx, ...
                    'EdgeColor', 'w', ...
                    'FaceColor', 'none', ...
                    'LineWidth', 1.5, ...
                    'LineStyle', '-');
                
                h2 = patch([acx acx acx acx], [acy acy acy acy], 'k',
                    'Parent', params.hAx, ...
                    'EdgeColor', params.mRbEdgeColor, ...
                    'FaceColor', params.mRbFaceColor, ...
                    'FaceAlpha', params.mRbFaceAlpha, ...
                    'LineWidth', 0.5, ...
                    'LineStyle', '-');    
 
                % create rubber band struct
                params.mRubberBand = struct(...
                    'obj',	[h1 h2], ...
                    'x1',  	acx, ...
                    'y1',  	acy, ...
                    'x2',  	acx, ...
                    'y2',  	acy);
                
                params.fIsRubberBandOn = true;
            end

        % Middle click: zooming
        case 'extend'
            %---ZoomMouseExtendBegin();
            if ~params.fIsZoomExtendAllowed
                
                % set new zoom grid for extend zoom
                [params.mZoomGrid, params.mZoomSteps] = ZoomLogGrid(params.mZoomMinPow, params.mZoomMaxPow, params.mZoomExtendNum);
                
                %---UpdateCurrentZoomAxes();
                [xLim, yLim] = GetAxesLimits(params.hAx);
                [curentZoomX, curentZoomY] = GetCurrentZoomAxesPercent(params.hAx, xLim, yLim, params.mDefaultXLim, params.mDefaultYLim);
                [nu, params.mZoomIndexX] = min(abs(params.mZoomGrid - curentZoomX));  %#ok ([~, ...])
                [nu, params.mZoomIndexY] = min(abs(params.mZoomGrid - curentZoomY));  %#ok ([~, ...])

                [wcx, wcy] = GetCursorCoordOnWindow(src);
                [acx, acy] = GetCursorCoordOnAxes(params.hAx);
                
                params.mStartX = wcx;
                params.mStartY = wcy;

                params.mBindX = acx;
                params.mBindY = acy;

                params.fIsZoomExtendAllowed = true;
            end
    end

    % Stores the updated parameter structure
    set(src, 'userdata', params);
end
%--------------------------------------------------------------------------

%==========================================================================
function WindowButtonUpCallback2D(src, evnt)      %#ok
    %WindowButtonUpCallback2D called when the mouse click is released.
    % Typically this is where we end the updating of the axis.
    
    params = get(src, 'userdata');

    if (isempty(params) || ~isfield(params, 'fIsDragAllowed'))
        return;
    end

    %%---DragMouseEnd();
    %if params.fIsDragAllowed
    %    params.fIsDragAllowed = false;
    %end

    %---ZoomMouseExtendEnd();
    if params.fIsZoomExtendAllowed
        %---SetDefaultZoomGrid();
        params.mZoomGrid = params.mDefaultZoomGrid;
        params.mZoomSteps = params.mDefaultZoomSteps;
        
        params.mZoomIndexX = find(params.mZoomGrid == 100);
        params.mZoomIndexY = params.mZoomIndexX;

        params.fIsZoomExtendAllowed = false;
    end

    %---RubberBandEnd();
    if params.fIsRubberBandOn
        params.fIsRubberBandOn = false;
        delete(params.mRubberBand.obj);          

        %---RubberBandZoomAxes();
        xLim = sort([params.mRubberBand.x1, params.mRubberBand.x2]);
        yLim = sort([params.mRubberBand.y1, params.mRubberBand.y2]);
        
        % Fix the final zoom factor to one of the values in the grid
        if (range(xLim) ~= 0 && range(yLim) ~= 0)
            [zoomPctX, zoomPctY] = GetCurrentZoomAxesPercent(params.hAx, ...
                xLim, yLim, params.mDefaultXLim, params.mDefaultYLim);
            
            if params.fIsImage
                zoomPctX = min(zoomPctX, zoomPctY);
                zoomPctY = zoomPctX;
            end
            
            cx = mean(xLim);
            cy = mean(yLim);
            
            xLim = RecalcZoomAxesLimits(xLim, params.mDefaultXLim, ...
                cx, zoomPctX, strcmp(get(params.hAx, 'xscale'), 'log'));
            yLim = RecalcZoomAxesLimits(yLim, params.mDefaultYLim, ...
                cy, zoomPctY, strcmp(get(params.hAx, 'yscale'), 'log'));
            
            SetAxesLimits(params.hAx, xLim, yLim);
        end

        params.mRubberBand = [];
    end

    % Store the updated parameters structure
    set(src, 'userdata', params);
end
%--------------------------------------------------------------------------

%==========================================================================
function WindowButtonMotionCallback2D(src, evnt)  %#ok
    %WindowButtonMotionCallback2D is called when the mouse is moved on
    % the figure. Typically this is where the axes are updated according
    % to the amplitude of the movement and update the position for the
    % next iteration.
            
    % Retrieve the parameter structure
    params = get(src, 'userdata');

    %---DragMouse();
    if params.fIsDragAllowed
        [cx, cy] = GetCursorCoordOnWindow(src);
        
        pdx = params.mStartX - cx;
        pdy = params.mStartY - cy;
        
        params.mStartX = cx;
        params.mStartY = cy;

        if (params.fIsImage)
            pdy = -pdy;
        end
        
        DragAxes(params.hAx, pdx, pdy, ...
            params.fIsEnableDragX, params.fIsEnableDragY);
    end

    %---RubberBandUpdate();
    if params.fIsRubberBandOn
        [acx, acy] = GetCursorCoordOnAxes(params.hAx);
        
        params.mRubberBand.x2 = acx;
        params.mRubberBand.y2 = acy;

        %---RubberBandSetPos();
        x1 = params.mRubberBand.x1;
        y1 = params.mRubberBand.y1;
        
        set(params.mRubberBand.obj, ...
            'XData', [x1 acx acx x1], ...
            'YData', [y1 y1 acy acy]);
    end
    
    %---ZoomMouseExtend();
    if params.fIsZoomExtendAllowed
        
        % Heuristic for pixel change to camera zoom factor 
        % (taken from function ZOOM, used in dragzoom.m)
        [wcx, wcy] = GetCursorCoordOnWindow(src);
        
        xy(1) = wcx - params.mStartX;
        xy(2) = wcy - params.mStartY;
        q = max(-0.9, min(0.9, sum(xy)/70)) + 1;

        % Move one step along the zooming grid
        if (q < 1)
            dz = -1;
        elseif (q > 1)
            dz = 1;
        else
            return;
        end
        
        %---ZoomAxes(direction, mZoom3DBindX, mZoom3DBindY)
        [xLim, yLim] = GetAxesLimits(params.hAx);
        
        % Keep a fixed axis ratio for images
        if params.fIsImage
            params.mZoomIndexX = params.mZoomIndexX + dz;
            params.mZoomIndexY = params.mZoomIndexY;
        else
            if params.fIsEnableZoomX
                params.mZoomIndexX = params.mZoomIndexX + dz;
            end
            if params.fIsEnableZoomY
                params.mZoomIndexY = params.mZoomIndexY + dz;
            end
        end

        % Make sure we stay within the zooming grid
        nz = length(params.mZoomGrid);
        params.mZoomIndexX = min(max(params.mZoomIndexX, 1), nz);
        params.mZoomIndexY = min(max(params.mZoomIndexY, 1), nz);

        xLim = RecalcZoomAxesLimits(xLim, params.mDefaultXLim, ...
            params.mBindX, params.mZoomGrid(params.mZoomIndexX), ...
            strcmp(get(params.hAx, 'xscale'), 'log'));
        yLim = RecalcZoomAxesLimits(yLim, params.mDefaultYLim, ...
            params.mBindY, params.mZoomGrid(params.mZoomIndexY), ...
            strcmp(get(params.hAx, 'yscale'), 'log'));
        
        SetAxesLimits(params.hAx, xLim, yLim);
        
        params.mStartX = wcx;
        params.mStartY = wcy;            
    end

    % Store the updated parameter structure
    set(src, 'userdata', params);
end
%--------------------------------------------------------------------------

%==========================================================================
function WindowKeyPressCallback2D(src, evnt)      %#ok
    %WindowKeyPressCallback2D is called when a keyboard key is pressed.
    % Typically used as an alterative to the mouse interactions and applies
    % discrete changes to the axis.

    % Get the parameter structure
    params = get(src, 'userdata');

    % We need a default value for the zoom factor to avoid duplicating
    % the rather length code for zooming.
    dz = 0;

    % Switch through the keys
    switch evnt.Key
        case {'0', 'equal'}
            %---ResetAxesToOrigView();
            SetAxesLimits(params.hAx, params.mDefaultXLim, params.mDefaultYLim);
            params.mZoomIndexX = find(params.mZoomGrid == 100);
            params.mZoomIndexY = params.mZoomIndexX;
        case {'add', 'plus', 'greater', 'comma'}
            %---ZoomKeys('plus');
            dz = 1;
        case {'hyphen', 'subtract', 'minus', 'less', 'period'}
            %---ZoomKeys('minus');
            dz = -1;
        case {'leftarrow', 'rightarrow', 'uparrow', 'downarrow', ...}
                'left', 'right', 'up', 'down'}
            %---DragKeys('...');
            dx = params.mDragShiftStep;
            dy = params.mDragShiftStep;
            
            % Increment of speed when you hold the button
            params.mDragShiftStep = params.mDragShiftStep + params.mDragShiftStepInc;
            
            % Determine which movement increment to keep
            switch evnt.Key(1)
                case 'r'
                    dy = 0;
                case 'l'
                    dx = -dx;
                    dy = 0;
                case 'd'
                    dy = -dy;
                    dx = 0;
                case 'u'
                    dx = 0;
            end

            % Y-axis is inverted in images
            if (params.fIsImage)
                dy = -dy;
            end

            DragAxes(params.hAx, dx, dy, ...
                params.fIsEnableDragX, params.fIsEnableDragY);
    
        case 'x'
            % Toggles zoom & drag restriction along the X axis
            % by forbidding movements along the Y axis
            params.fIsEnableDragY = ~params.fIsEnableDragY;
            params.fIsEnableZoomY = ~params.fIsEnableZoomY;

        case 'y'
            % Toggles along the Y axis
            params.fIsEnableDragX = ~params.fIsEnableDragX;
            params.fIsEnableZoomX = ~params.fIsEnableZoomX;
    end

    % Merged both zoom actions into a single one
    if (dz ~= 0)
        %---UpdateCurrentZoomAxes();
        
        % set new zoom grid for extend zoom
        [params.mZoomGrid, params.mZoomSteps] = ZoomLogGrid(params.mZoomMinPow, params.mZoomMaxPow, params.mZoomExtendNum);
        
        %---UpdateCurrentZoomAxes();
        [xLim, yLim] = GetAxesLimits(params.hAx);
        [curentZoomX, curentZoomY] = GetCurrentZoomAxesPercent(params.hAx, xLim, yLim, params.mDefaultXLim, params.mDefaultYLim);
        [nu, params.mZoomIndexX] = min(abs(params.mZoomGrid - curentZoomX));  %#ok ([~, ...])
        [nu, params.mZoomIndexY] = min(abs(params.mZoomGrid - curentZoomY));  %#ok ([~, ...])

        [acx, acy] = GetCursorCoordOnAxes(params.hAx);
        
        if params.fIsImage
            params.mZoomIndexX = params.mZoomIndexX + dz;
            params.mZoomIndexY = params.mZoomIndexX;
        else
            if params.fIsEnableZoomX
                params.mZoomIndexX = params.mZoomIndexX + dz;
            end
            if params.fIsEnableZoomY
                params.mZoomIndexY = params.mZoomIndexY + dz;
            end
        end

        % Make sure we stay within the zooming grid
        nz = length(params.mZoomGrid);
        params.mZoomIndexX = min(max(params.mZoomIndexX, 1), nz);
        params.mZoomIndexY = min(max(params.mZoomIndexY, 1), nz);

        xLim = RecalcZoomAxesLimits(xLim, params.mDefaultXLim, ...
            acx, params.mZoomGrid(params.mZoomIndexX), ...
            strcmp(get(params.hAx, 'xscale'), 'log'));
        yLim = RecalcZoomAxesLimits(yLim, params.mDefaultYLim, ...
            acy, params.mZoomGrid(params.mZoomIndexY), ...
            strcmp(get(params.hAx, 'yscale'), 'log'));
        
        SetAxesLimits(params.hAx, xLim, yLim);
    end

    % Store the updated parameter structure
    set(src, 'userdata', params);
    drawnow;
end
%--------------------------------------------------------------------------

%==========================================================================
function WindowKeyReleaseCallback2D(src, evnt)    %#ok
    %WindowKeyReleaseCallback2D is called when a pressed key is released.
    % The only usage is to restore the incremented drag shift to its
    % original value.

    switch evnt.Key
        case {'leftarrow', 'rightarrow', 'uparrow', 'downarrow', ...}
                'left', 'right', 'up', 'down'}
            params = get(src, 'userdata');
            params.mDragShiftStep = params.mDragSaveShiftStep;
            set(src, 'userdata', params);
    end
end
%--------------------------------------------------------------------------

%==========================================================================
function DragAxes(hAx, pdx, pdy, fIsEnableDragX, fIsEnableDragY)
    %DragAxes a subfunction responsible for calculating the new axes limits
    % based on the provided movement increment.
    
    [xLim, yLim] = GetAxesLimits(hAx);
    
    %---pos = GetObjPos(hAx, 'Pixels');
    dfltUnits = get(hAx, 'Units');
    set(hAx, 'Units', 'Pixels');
    pos = get(hAx, 'Position');
    set(hAx, 'Units', dfltUnits);

    pbar = get(hAx, 'PlotBoxAspectRatio');
    
    % Heuristic for replacing the PAN function that has some type
    % of bug according to Evgeny (used in dragzoom.m)
    imAspectRatioX = pbar(2) / pbar(1);
    if (imAspectRatioX ~= 1)
        posAspectRatioX = pos(3) / pos(4);
        arFactorX = imAspectRatioX * posAspectRatioX;
        if (arFactorX < 1)
            arFactorX = 1;
        end
    else
        arFactorX = 1;
    end
    
    imAspectRatioY = pbar(1) / pbar(2);
    if (imAspectRatioY ~= 1)
        posAspectRatioY = pos(4) / pos(3);
        arFactorY = imAspectRatioY * posAspectRatioY;
        if (arFactorY < 1)
            arFactorY = 1;
        end
    else
        arFactorY = 1;
    end
    
    if fIsEnableDragX
        % For log plots, transform to linear scale
        if strcmp(get(hAx, 'xscale'), 'log')
            xLim = log10(xLim);
            xLim = FixInfLogLimits('x', xLim);
            isXLog = true;
        else
            isXLog = false;
        end
        
        dx = pdx * range(xLim) / (pos(3) / arFactorX);
        xLim = xLim + dx;
        
        % For log plots, untransform limits
        if isXLog
            xLim = 10.^(xLim);
        end
    end
    if fIsEnableDragY
        if strcmp(get(hAx, 'yscale'), 'log')
            yLim = log10(yLim);
            yLim = FixInfLogLimits('y', yLim);
            isYLog = true;
        else
            isYLog = false;
        end
        
        dy = pdy * range(yLim) / (pos(4) / arFactorY);
        
        yLim = yLim + dy; 
        
        if isYLog
            yLim = 10.^(yLim);
        end
    end
    
    % Applies the new limits to the axes
    SetAxesLimits(hAx, xLim, yLim);
end
%--------------------------------------------------------------------------

%==========================================================================
function [curentZoomX, curentZoomY] = GetCurrentZoomAxesPercent(hAx, xLim, yLim, mDefaultXLim, mDefaultYLim)
    %GetCurrentZoomAxesPercent determines given the current axes limits 
    % which zoom factor we are using.
    
    if strcmp(get(hAx, 'xscale'), 'log')
        xLim = log10(xLim);
        defaultXLim = log10(mDefaultXLim);
    else
        defaultXLim = mDefaultXLim;
    end
    if strcmp(get(hAx, 'yscale'), 'log')
        yLim = log10(yLim);
        defaultYLim = log10(mDefaultYLim);
    else
        defaultYLim = mDefaultYLim;
    end
    
    curentZoomX = range(defaultXLim) * 100 / range(xLim);
    curentZoomY = range(defaultYLim) * 100 / range(yLim);
end
%--------------------------------------------------------------------------

%==========================================================================
function [x, y] = GetCursorCoordOnAxes(hAx)
    %GetCursorCoordOnAxImg helper function to get the position of the
    % mouse in the coordinates of the axis.
    
    crd = get(hAx, 'CurrentPoint');
    x = crd(2,1);
    y = crd(2,2);
end
%--------------------------------------------------------------------------

%==========================================================================
function [x, y] = GetCursorCoordOnWindow(hFig)
    %GetCursorCoordOnWindow get the position of the mouse on the figure
    % in pixels
    
    %dfltUnits = get(hFig, 'Units');
    %set(hFig, 'Units', 'pixels');
    
    crd = get(hFig, 'CurrentPoint');
    x = crd(1); 
    y = crd(2);
    
    %set(hFig, 'Units', dfltUnits);
end
%--------------------------------------------------------------------------

%==========================================================================
function [xLim, yLim] = GetAxesLimits(hAx)
    %GetAxesLimits
    
    xLim = get(hAx, 'XLim');
    yLim = get(hAx, 'YLim');
end
%--------------------------------------------------------------------------

%==========================================================================
function SetAxesLimits(hAx, xLim, yLim)
    %SetAxesLimits
    
    set(hAx, 'XLim', xLim);
    set(hAx, 'YLim', yLim);
end
%--------------------------------------------------------------------------

%==========================================================================
function axLim = RecalcZoomAxesLimits(axLim, axLimDflt, zcCrd, zoomPct, isLog)
    %RecalcZoomAxesLimits recalculates the axes limits
    
    if isLog
        axLim = log10(axLim);
        %---axLim = FixInfLogLimits(ax, axLim);
        axLimDflt = log10(axLimDflt);
        zcCrd = log10(zcCrd);

        % Simple hack to avoid using the original undocumented
        % axis function
        if (~all(isfinite(axLim)) || ~all(isreal(axLim)))
            axLim = axLimDflt;
        end
    end
            
    if (zcCrd < axLim(1)), zcCrd = axLim(1); end
    if (zcCrd > axLim(2)), zcCrd = axLim(2); end
    
    rf = range(axLim);
    ra = range([axLim(1), zcCrd]);
    rb = range([zcCrd, axLim(2)]);
    
    cfa = ra / rf; 
    cfb = rb / rf;
    
    newRange = range(axLimDflt) * 100 / zoomPct;
    dRange = newRange - rf;
    
    axLim(1) = axLim(1) - dRange * cfa;
    axLim(2) = axLim(2) + dRange * cfb;
    
    if isLog
        axLim = 10.^axLim;
    end
end
%--------------------------------------------------------------------------

%==========================================================================
function [zg, st] = ZoomLogGrid(a, b, n)
    %ZoomLogGrid creates the log zoom grid
    
    zg = unique(round(logspace(a, b, n)));
    
    zg(zg<10) = [];	% begin zoom == 10%
    st = length(zg);
    
end
%--------------------------------------------------------------------------
