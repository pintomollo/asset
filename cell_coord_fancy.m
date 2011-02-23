function cell_coord_fancy(mymovie)

close all hidden;

ierr = 0;
%maxint = 1;
hwait = waitbar(0,' ','Name','CellCoord Info','Tag','wait','Visible','off');
set(hwait,'HandleVisibility','on');

% Change the root default font size for question dialog boxes.

set(0,'DefaultTextFontSize',14)
%hwait = waitbar(0,' ','Name','CellCoord Info','Tag','wait','Visible','off');

if(~isfield(mymovie,'edges') && (~isfield(mymovie,'egg') || ~isfield(mymovie,'cortex')))
  %mymovie = compute_edges(mymovie);

  %movefile(mymovie.edges,['edges' mymovie.edges]);
  %mymovie.edges = ['edges' mymovie.edges];

  %mymovie = extract_pts(mymovie,10);

  %save 'GZ869-sp5.mat' mymovie;
end
if(~isfield(mymovie,'empty'))
  mymovie.empty = 'empty.tmp';

  %save 'GZ869-sp5.mat' mymovie;
end

[imgsize, nframes] = size_data(mymovie.dic);

%mov_data = struct([]);
hwindows = struct([]);
%user_interf = struct([]);

%imgsize = size(mymovie(1).data);
%nframes = size(mymovie,2);
mov_data = struct('min_dist',max(imgsize) / 15, 'nframes', nframes, 'nsplines', 0, 'splines', [], 'ref_spline', []);
%mov_data.min_dist = max(imgsize) / 15;
%mov_data.nframes = nframes;
%mov_data.nsplines = 0;
%mov_data.splines = [];

user_interf = struct('current_frame', 1, 'previous_frame', 1, 'selected', 0, 'drag_cp', false, 'operation', '', 'index_main', 0, 'index_ref', 0, 'upmax', 20, 'ms2',12, 'sl_color', [0 0 1], 'old_frac', 0, 'old_lvl', 0);
%user_interf.current_frame = 1; 
%user_interf.previous_frame = 1;
%user_interf.selected = 0;
%user_interf.drag_cp = false;
%user_interf.operation = '';

% Create the primary figure and axes.

hwindows = create_window('main',mymovie.dic, 'DIC');
user_interf.index_main = length(hwindows);
user_interf.old_frac = get(findobj('Tag','uifrac'),'Value');
user_interf.old_lvl = get(findobj('Tag','uidetectlvl'),'Value');
hwindows = create_window('disp',mymovie.c010, 'PAR-6');
%if(~isfield(mymovie,'c001'))
%  hwindows = create_window('disp',mymovie.c100, 'PAR-6');
%end

hwindows = create_window('coord',mymovie.empty, '');
user_interf.index_ref = length(hwindows);
%hwindows = create_window('evol',mymovie.edges);

htimer = timer('TimerFcn', @NextFrame, ...
          'Period', 0.25, ...
          'ExecutionMode', 'fixedDelay', ...
          'BusyMode', 'queue');

create('egg');
%tscurve;

clear ierr imgsize nframes;

return;


%***********************************************************
% TSPACKGUI Nested Functions
%***********************************************************

function create(obj)
% Create line objects for a new data set.
%
% USAGE:  createlo
%
% Adds line objects harrow, hcp0, hcpi0, hcpoly, hcs0, and 
% hcsp to the axes in figure 1, adds hcp1 and hcs1 to 
% figure 2, and adds hcp2 and hcs2 to figure 3.  In figure 
% 1, each quiver3 object, curve segment, and control point 
% is an individually selectable line object.

  hax = hwindows(user_interf.index_main).hax;
  mov_data.nsplines = mov_data.nsplines + 1;

  hcp = [];
  hcs = [];

  frac = get(findobj('Tag','uifrac'),'Value');

switch obj
case 'egg'
  cp_color = [0 1 0];
  cs_color = [0 1 0];

  if(user_interf.index_ref~=0)
    tmpx = get(hwindows(user_interf.index_ref).hax, 'XLim');
    tmpy = get(hwindows(user_interf.index_ref).hax, 'YLim');

    eggsize = [tmpx(2) tmpy(2)] / 1.1;
    %eggsize = [50 30] / 2;

    [eggx, eggy] =  draw_ellipse([0,0],eggsize,0,8);

    axex = [1 0 -1 0] * 1.1;
    axey = [0 1 0 -1] * 1.1;

    %eggx = eggx * eggsize(1);
    %eggy = eggy * eggsize(2);

    axex = axex * eggsize(1);
    axey = axey * eggsize(2);

    %eggx = eggx + eggsize(1) + 5;
    %eggy = eggy + eggsize(2) + 5;

    %axex = axex + eggsize(1) + 5;
    %axey = axey + eggsize(2) + 5;

    %t = arcl2d(eggx,eggy);

    %egg = mycsape(t, [eggx;eggy], 'periodic');

    mov_data.ref_spline.x = eggx(1:end);
    mov_data.ref_spline.y = eggy(1:end);
   
   egghcs = line(eggx(1),eggy(1), ...
                  'Color',cs_color, ...
                  'HandleVisibility', 'callback', ...
                  'EraseMode', 'none', ...
                  'HitTest','off', ...
                  'LineStyle','-', 'Marker','none', ...
                  'Parent',hwindows(user_interf.index_ref).hax, ...
                  'SelectionHighlight','off', ...
                  'Tag',num2str(mov_data.nsplines));
    hwindows(user_interf.index_ref).hsplines(mov_data.nsplines).hcs = egghcs;

    [h1,h2,hlegends,legends] = legend(hwindows(user_interf.index_ref).hax);
    legends = cat(1,legends,'Eggshell');
    hlegends = cat(1,hlegends,egghcs);
    legend(hwindows(user_interf.index_ref).hax, hlegends, legends, 'Location', [0.05 0.7 0.1 0.1] );

    legend(hwindows(user_interf.index_ref).hax,'boxoff');

    line(axex([1 3]), axey([1 1]), ...
            'Color', 'k', ...
            'HandleVisibility', 'callback', ...
            'EraseMode', 'none', ...
            'HitTest','off', ...
            'LineStyle','-', 'Marker','none', ...
            'Parent',hwindows(user_interf.index_ref).hax, ...
            'SelectionHighlight','off');

    line(axex([2 2]), axey([2 4]), ...
            'Color', 'k', ...
            'HandleVisibility', 'callback', ...
            'EraseMode', 'none', ...
            'HitTest','off', ...
            'LineStyle','-', 'Marker','none', ...
            'Parent',hwindows(user_interf.index_ref).hax, ...
            'SelectionHighlight','off');

      text(axex + [2 0 -2 0], axey + [0 2 0 -2], {'0'; '3\pi/2'; '\pi'; '\pi/2'});
  end

 %mov_data.ref_spline.spline = egg;

  %img = load_data(mymovie.edges,user_interf.current_frame);
  %[x,y] = detect_border(img,[],get(findobj('Tag','uidetectlvl'),'Value'),true);
  %[x,y,t] = optimize_spline(x,y,mov_data.min_dist);
  %[x,y,t] = optimize_spline(x,y,img,frac,0.01);
  x = mymovie.egg(user_interf.current_frame).path(:,1);
  y = mymovie.egg(user_interf.current_frame).path(:,2);
  n = length(x);
  x = x(:)';
  y = y(:)';
  hcp = zeros(1,n);
  hcs = zeros(1,n);

  for i = 1:n

    % Create control point line objects.

     hcp(i) = line(x(i),y(i), ...
                  'ButtonDownFcn',@CpButtonDnFcn, ...
                  'HandleVisibility','callback', ...
                  'HitTest','off', ...
                  'Visible','off', ...
                  'LineStyle','none','Marker','o', ...
                  'MarkerEdgeColor',cp_color, ...
                  'MarkerFaceColor',cp_color, ...
                  'MarkerSize',user_interf.ms2, ...
                  'Parent',hax, ...
                  'SelectionHighlight','off', ...
                  'Tag',num2str(i));
  end
   hcs = line(x(1),y(1), ...
                  'ButtonDownFcn',@CsButtonDnFcn, ...
                  'Color',cs_color, ...
                  'HandleVisibility', 'callback', ...
                  'EraseMode', 'none', ...
                  'HitTest','off', ...
                  'LineStyle','-', 'Marker','none', ...
                  'Parent',hax, ...
                  'SelectionHighlight','off', ...
                  'Tag',num2str(mov_data.nsplines));
  case 'cortex'

        %img = get(himg,'CData');
        %if(size(img,3)>1)
        %  img = rgb2gray(img);
        %end
    %img = load_data(mymovie.edges,user_interf.current_frame);
    %egg_mask = create_mask(mov_data.splines(1,user_interf.current_frame).ppform,img);
    %egg_mask = imerode(egg_mask,strel('disk',3));
    %img = img.*egg_mask;

    %[x,y] = detect_border(img,mov_data.splines(1,user_interf.current_frame).ppform,get(findobj('Tag','uidetectlvl'),'Value'),true);
    %[x,y,t] = optimize_spline(x,y,img,frac,0.01);

    x = mymovie.cortex(user_interf.current_frame).path(:,1);
    y = mymovie.cortex(user_interf.current_frame).path(:,2);

    % Convert column vectors to row vectors.
    n = length(x);
    x = x(:)';
    y = y(:)';
    %hcp = zeros(1,n-1);
    %hcs = zeros(1,n-1);
    cp_color = [1 0.66 0];
    cs_color = [1 0.66 0];

        for i=1:n
         hcp(1) = line(x(1),y(1), ...
                        'ButtonDownFcn',@CpButtonDnFcn, ...
                        'HandleVisibility','callback', ...
                        'HitTest','off', ...
                        'Visible','off', ...
                        'LineStyle','none','Marker','o', ...
                        'MarkerEdgeColor',cp_color, ...
                        'MarkerFaceColor',cp_color, ...
                        'MarkerSize',user_interf.ms2, ...
                        'Parent',hax(1), ...
                        'SelectionHighlight','off', ...
                    'Tag',num2str(i));
        end


      hcs = line(x(1),y(1), ...
                  'ButtonDownFcn',@CsButtonDnFcn, ...
                  'Color',cs_color, ...
                  'HandleVisibility', 'callback', ...
                  'HitTest','off', ...
                  'LineStyle','-', 'Marker','none', ...
                  'Parent',hax, ...
                  'SelectionHighlight','off', ...
                  'Tag',num2str(mov_data.nsplines)); 

          if(user_interf.index_ref~=0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%% USE SURF INSTEAD OF SPLINES %%%%%%%%%%%%%%

           eggx = mov_data.ref_spline.x(1);
           eggy = mov_data.ref_spline.y(1);
           
           egghcs = line(eggx,eggy, ...
                          'Color',cs_color, ...
                          'HandleVisibility', 'callback', ...
                          'EraseMode', 'none', ...
                          'HitTest','off', ...
                          'LineStyle','-', 'Marker','none', ...
                          'Parent',hwindows(user_interf.index_ref).hax, ...
                          'SelectionHighlight','off', ...
                          'Tag',num2str(mov_data.nsplines));

            hwindows(user_interf.index_ref).hsplines(mov_data.nsplines).hcs = egghcs;

            [h1,h2,hlegends,legends] = legend(hwindows(user_interf.index_ref).hax);
            legends = cat(1,legends,'Cortex');
            hlegends = cat(1,hlegends,egghcs);
            legend(hwindows(user_interf.index_ref).hax, hlegends, legends, 'Location', [0.05 0.7 0.1 0.1] );
            legend(hwindows(user_interf.index_ref).hax,'boxoff');
          end
  end


  new_splines = repmat(struct('cp_color',[],'cs_color',[], ...
                   'ppform',[],'x',[],'y',[], 'dist', [], 'profile', [], ...
                   'ulist',[],'ulcnt',0,'uptr',0, 'name',''),mov_data.nsplines,mov_data.nframes);
  new_hsplines = repmat(struct('hcp',[],'hcs',[], ...
                   'n',0,'t',[]),1,mov_data.nsplines);

  if(mov_data.nsplines > 1)
    new_splines(1:mov_data.nsplines-1,:) = mov_data.splines;
    new_hsplines(1:mov_data.nsplines-1) = hwindows(user_interf.index_main).hsplines(:);
  end

  mov_data.splines = new_splines;
  hwindows(user_interf.index_main).hsplines = new_hsplines;

  mov_data.splines(mov_data.nsplines,user_interf.current_frame).name = obj;
  mov_data.splines(mov_data.nsplines,user_interf.current_frame).x = x(:)';
  mov_data.splines(mov_data.nsplines,user_interf.current_frame).y = y(:)';
  mov_data.splines(mov_data.nsplines,user_interf.current_frame).ulist = repmat(struct('opcode',' ', ...
                                    'index',0, 'data',[]), 1, user_interf.upmax);
  mov_data.splines(mov_data.nsplines,user_interf.current_frame).cp_color = cp_color;
  mov_data.splines(mov_data.nsplines,user_interf.current_frame).cs_color = cs_color;

  %hsplines(nsplines).t = t;
  %hsplines(nsplines).n = n;
  %hsplines(nsplines).hcp = hcp;
  hwindows(user_interf.index_main).hsplines(mov_data.nsplines).hcs = hcs;
  hwindows(user_interf.index_main).hsplines(mov_data.nsplines).hcp = hcp;

  user_interf.selected = mov_data.nsplines;
  tsknots;
  tscurve;
  user_interf.selected = 0;
           
  return;
end  % createlo

%***********************************************************
function [new_hwindows] = create_window(type, img_path, fig_name)

  hwindow = struct('hsplines', [], 'hc', [], 'hdragcp', [], 'hfig', [], 'hui', [], 'hgroup', [], 'hax', [], 'himg', [], 'img_path', img_path);
  new_hwindows = repmat(hwindow,1,length(hwindows)+1);

  hfig = figure;
  drawnow  % This may be required to correctly set properties.
  ssiz = get(0,'ScreenSize');

  hui = [];
  hgroup = [];
  htmp = [];

  pos = 10*ones(2,2);
  spacing = 2;
  rowdist = 20;
  for i=1:size(pos,2)
    pos(i,2) = spacing*i + rowdist*(i-1);
  end

  axes_position = [0 0.05 1 0.95];

  [imgsize, nframes] = size_data(img_path);
  ylimsize = [0 imgsize(1)];
  xlimsize = [0 imgsize(2)];
  
  switch type
    case 'main'

      if(nargin<3)
        fig_name = 'CellCoord';
      end

      if(nframes==1)
        slidernum = 1;
      else
        slidernum = nframes-1;
      end

      set(hfig, 'CloseRequestFcn',@FexitFcn, ...
          'DoubleBuffer','on', ...
          'KeyPressFcn',@KeyPressFcn, ...
          'MenuBar','none', ...
          'Name',fig_name, 'NumberTitle','off', ...
          'Units','pixels', ...
          'Position',[10 ssiz(4)/2 ssiz(3)/2 ssiz(4)/2], ...
          'Toolbar','none', ...
          'Tag', type, ...
          'WindowButtonMotionFcn',@MotionFcn, ...
          'WindowButtonUpFcn',@ButtonUpFcn)
      bcolor = get(hfig,'Color');

      % Create a group of handlers
      hgroup(end+1).tag = 'video';
      hgroup(end+1).hui = [];

      uisize = [181 22];
      rowindx = 2;
      htmp(end+1) = uicontrol('Parent',hfig, ...
            'Position',[pos(rowindx,:) uisize], ...
          'Callback', @SelectFrame, ...
          'Min',1, 'Max',nframes, ...
          'SliderStep',[1/slidernum 5/slidernum], ...
          'Style','slider', ...
          'Tag','uiframes', ...
          'Value', 1 );
      pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      uisize = [45 15];
      rowindx = 2;
      htmp(end+1) = uicontrol('Parent', hfig, ...
            'Position',[pos(rowindx,:) uisize], ...
          'Tag','uiframestext', ...
          'Style', 'text', ...
          'FontWeight','Bold', ...
          'BackgroundColor', bcolor, ...
          'String', sprintf('%i / %i', user_interf.current_frame, nframes));
      pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      uisize = [50 20];
      rowindx = 1;
      htmp(end+1) = uicontrol('Parent',hfig, ...
            'Position',[pos(rowindx,:) uisize], ...
          'Callback', @PlayMovie, ...
          'FontWeight','Bold', ...
          'String','Play', ...
          'Tag', 'uiplay', ...
          'Style','togglebutton');
      pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      hgroup(end).hui = htmp;
      htmp = [];

      uisize = [80 20];
      rowindx = 1;
      hui(end+1) = uicontrol('Parent',hfig, ...
            'Position',[pos(rowindx,:) uisize], ...
            'FontWeight','Bold', ...
            'Callback', @AdaptSpline, ...
            'Value', 0, ...
            'String','Adapt Spline', ...
            'Style','togglebutton', ...
            'Tag', 'uiadapt');
      pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      uisize = [80 20];
      rowindx = 1;
      hui(end+1) = uicontrol('Parent',hfig, ...
            'Position',[pos(rowindx,:) uisize], ...
            'Callback', @Redetect, ...
            'FontWeight','Bold', ...
            'String','Redetect', ...
            'Style','pushbutton', ...
            'Tag', 'uiprofile');
      pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      %uisize = [60 15];
      %rowindx = 1;
      %htmp(end+1) = uicontrol('Parent', hfig, ...
      %      'Position',[pos(rowindx,:) uisize], ...
      %    'Style', 'text', ... 
      %    'FontWeight','Bold', ...
      %    'BackgroundColor', bcolor, ...
      %    'String', 'Copy Prev.');
      %pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      %uisize = [20 20];
      %rowindx = 1;
      %htmp(end+1) = uicontrol('Parent',hfig, ...
      %      'Position',[pos(rowindx,:) uisize], ...
      %    'BackgroundColor', bcolor, ...
      %    'Style','checkbox', ...
      %    'Value', true, ...
      %    'Tag', 'uicopy');
      %pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      uisize = [35 15];
      rowindx = 1;
       hui(end+1) = uicontrol('Parent', hfig, ...
            'Position',[pos(rowindx,:) uisize], ...
          'Style', 'text', ... 
          'FontWeight','Bold', ...
          'BackgroundColor', bcolor, ...
          'String', 'Cortex');
      pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      uisize = [20 20];
      rowindx = 1;
      hui(end+1) =  uicontrol('Parent',hfig, ...
            'Position',[pos(rowindx,:) uisize], ...
          'Callback',@DetectElement, ...
          'BackgroundColor', bcolor, ...
          'Style','checkbox', ...
          'Tag', 'uicortex');
      pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      uisize = [35 15];
      rowindx = 1;
       hui(end+1) = uicontrol('Parent', hfig, ...
            'Position',[pos(rowindx,:) uisize], ...
          'Style', 'text', ... 
          'FontWeight','Bold', ...
          'BackgroundColor', bcolor, ...
          'String', 'A<->P');
      pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      uisize = [20 20];
      rowindx = 1;
      hui(end+1) =  uicontrol('Parent',hfig, ...
            'Position',[pos(rowindx,:) uisize], ...
          'Callback',@InvertAP, ...
          'BackgroundColor', bcolor, ...
          'Style','checkbox', ...
          'Tag', 'uiinvert');
      pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;


      %uisize = [80 15];
      %rowindx = 2;
      % hui(end+1) =uicontrol('Parent', hfig, ...
      %      'Position',[pos(rowindx,:) uisize], ...
      %    'Style', 'text', ...
      %    'FontWeight','Bold', ...
      %    'BackgroundColor', bcolor, ...
      %    'String', 'Detection Level');
      %pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      %uisize = [181 22];
      %rowindx = 2;
      %hui(end+1) =uicontrol('Parent',hfig, ...
      %      'Position',[pos(rowindx,:) uisize], ...
      %      'Callback', @Redetect, ...
      %      'Min',0, 'Max', 2, ...
      %      'SliderStep',[1/100 1/100], ...
      %      'Style','slider', ...
      %      'Tag','uidetectlvl', ...
      %      'Value', 1 );
      %pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      %uisize = [80 15];
      %rowindx = 2;
      % hui(end+1) =uicontrol('Parent', hfig, ...
      %      'Position',[pos(rowindx,:) uisize], ...
      %    'Style', 'text', ...
      %    'FontWeight','Bold', ...
      %    'BackgroundColor', bcolor, ...
      %    'String', 'Roundness Infl.');
      %pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      %uisize = [90 22];
      %rowindx = 2;
      %hui(end+1) =uicontrol('Parent',hfig, ...
      %      'Position',[pos(rowindx,:) uisize], ...
      %      'Callback', @Redetect, ...
      %      'Min',0, 'Max', 1, ...
      %      'SliderStep',[1/20 1/20], ...
      %      'Style','slider', ...
      %      'Tag','uifrac', ...
      %      'Value', 0.6 );
      %pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      %uisize = [80 15];
      %rowindx = 1;
      % hui(end+1) =uicontrol('Parent', hfig, ...
      %      'Position',[pos(rowindx,:) uisize], ...
      %    'Style', 'text', ...
      %    'FontWeight','Bold', ...
      %    'BackgroundColor', bcolor, ...
      %    'String', 'Init. Temp');
      %pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      %uisize = [90 22];
      %rowindx = 1;
      %hui(end+1) =uicontrol('Parent',hfig, ...
      %      'Position',[pos(rowindx,:) uisize], ...
      %      'Callback', @Redetect, ...
      %      'Style','edit', ...
      %      'Tag','uitemp', ...
      %      'String', num2str(0.75) );
      %pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      if(nframes == 1)
        set(hgroup(end).hui, 'Visible', 'off');
      end

      for i=1:length(hgroup)
        hui = [hui hgroup(i).hui];
      end
      %for i=1:length(hgroup)
      %  hgroup(i).hnot = hui(~hgroup(i).hui);
      %end

      % File Menu

      hmenu = uimenu('Label','File', ...
                    'Tag','uimenufile', ...
                     'Parent',hfig);
       uimenu('Accelerator','o', ...
                      'Callback',@FopenFcn, ...
                      'Label','Open', ...
                      'Parent',hmenu);
       uimenu('Accelerator','s', ...
                      'Callback',@FsaveFcn, ...
                      'Label','Save', ...
                      'Parent',hmenu);
       uimenu('Callback',@FexitFcn, ...
                      'Label','Exit', ...
                      'Parent',hmenu);
      % Operations Menu

      hmenu = uimenu('Label','Operations', ...
                     'Tag', 'uimenuop', ...
                     'Enable','off', ...
                     'Parent',hfig);
       uimenu('Accelerator','v', ...
              'Callback',@MoveCpFcn, ...
              'Label','Move Control Point', ...
              'Parent',hmenu)
       uimenu('Accelerator','i', ...
              'Callback',@InsertCpFcn, ...
              'Label','Insert Control Point', ...
              'Parent',hmenu)
       uimenu('Accelerator','d', ...
              'Callback',@DeleteCpFcn, ...
              'Label','Delete Control Point', ...
              'Parent',hmenu)
       uimenu('Accelerator','u', ...
              'Callback',@UndoFcn, ...
              'Label','Undo Last Operation', ...
              'Parent',hmenu, 'Separator','on')
       uimenu('Accelerator','s', ...
              'Callback',@SelectFcn, ...
              'Label','Select a spline', ...
              'Parent',hmenu)

    case 'disp'

      if(nargin<3)
        fig_name = 'CellCoord Channel';
      end

      nlayers = findobj('Tag',type);
      nlayers = length(nlayers);

      set(hfig, 'CloseRequestFcn',@FvoidFcn, ...
          'DoubleBuffer','on', ...
          'MenuBar','none', ...
          'Name',fig_name, 'NumberTitle','off', ...
          'Units','pixels', ...
          'Tag', type, ...
          'Position',[(100*(nlayers+1)) + ssiz(4)/2 ssiz(4)/2 ssiz(3)/2 ssiz(4)/2], ...
          'Toolbar','none');
      bcolor = get(hfig,'Color');

      %uisize = [80 20];
      %rowindx = 2;
      %hui(end+1) = uicontrol('Parent',hfig, ...
      %      'Position',[pos(rowindx,:) uisize], ...
      %      'Callback', @FitSpline, ...
      %      'FontWeight','Bold', ...
      %      'String','Fit', ...
      %      'Style','pushbutton', ...
      %      'Tag', 'uifit');
      %pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      %uisize = [80 20];
      %rowindx = 1;
      %hui(end+1) = uicontrol('Parent',hfig, ...
      %      'Position',[pos(rowindx,:) uisize], ...
      %      'Callback', @Score, ...
      %      'FontWeight','Bold', ...
      %      'String','Score', ...
      %      'Style','pushbutton', ...
      %      'Tag', 'uiscore');
      %pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      %uisize = [150 15];
      %rowindx = 1;
      %hui(end+1) = uicontrol('Parent', hfig, ...
      %      'Position',[pos(rowindx,:) uisize], ...
      %    'Style', 'text', ... 
      %    'FontWeight','Bold', ...
      %    'BackgroundColor', bcolor, ...
      %    'Tag','uiscoretxt',...
      %    'String', '');
      %pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      %uisize = [80 20];
      %rowindx = 2;
      %hui(end+1) = uicontrol('Parent',hfig, ...
      %      'Position',[pos(rowindx,:) uisize], ...
      %      'FontWeight','Bold', ...
      %      'String',{'radial', 'neigborhood', 'montecarlo'}, ...
      %      'Style','popupmenu', ...
      %      'Tag', 'uifittype');
      %pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

  case 'evol'

      set(hfig, 'CloseRequestFcn',@FvoidFcn, ...
          'DoubleBuffer','on', ...
          'MenuBar','none', ...
          'Name','Simulated Annealing', 'NumberTitle','off', ...
          'Units','pixels', ...
          'Tag', type, ...
          'Position',[ssiz(4)/2 ssiz(4)/2-200 ssiz(3)/2 ssiz(4)/2], ...
          'Toolbar','none');
      bcolor = get(hfig,'Color');

      uisize = [80 15];
      rowindx = 2;
       hui(end+1) =uicontrol('Parent', hfig, ...
            'Position',[pos(rowindx,:) uisize], ...
          'Style', 'text', ...
          'Visible', 'off', ...
          'Tag', 'filename', ...
          'String', '');
      pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

  case 'coord'

      set(hfig, 'CloseRequestFcn',@FvoidFcn, ...
          'DoubleBuffer','on', ...
          'MenuBar','none', ...
          'Name','Reference System', 'NumberTitle','off', ...
          'Units','pixels', ...
          'Tag', type, ...
          'Position',[ssiz(4)/2 ssiz(4)/2-200 ssiz(3)/1.5 ssiz(4)/2], ...
          'Toolbar','none');
      bcolor = get(hfig,'Color');

      uisize = [80 20];
      rowindx = 1;
      hui(end+1) = uicontrol('Parent',hfig, ...
            'Position',[pos(rowindx,:) uisize], ...
            'Callback', @GetProfiles, ...
            'FontWeight','Bold', ...
            'String','Get Profiles', ...
            'Style','pushbutton', ...
            'Tag', 'uiprofile');
      pos(rowindx,1) = pos(rowindx,1) + uisize(1) + spacing;

      imgsize = [40 60];

      hax2 = axes;
      set(hax2, ...
        'FontSize',18, ...
        'XLim', [0 2*pi], ...
        'XLimMode', 'manual', ...
        'XTick', [0:pi/2:2*pi], ...
        'XTickLabel', {'0'; 'pi/2'; 'pi'; '3*pi/2'; '2*pi'}, ...
        'YLim', [0 1], ...
        'YLimMode', 'manual', ...
        'YTick', [0:0.1:1], ...
        'NextPlot','add', ...
        'OuterPosition', [0.5 0.05 0.5 0.95], ...
        'Tag', 'profile_plot',...
        'Parent',hfig);
        %'YTickLabel', {'0'; '10'; '20'; '30'; '40'; '50'; '60'; '70'; '80'; '90'; '100'}, ...
      %'XTick', [0:0.1:1], ...
      %'XTickLabel', {'0'; '10'; '20'; '30'; '40'; '50'; '60'; '70'; '80'; '90'; '100'}, ...

      xlabel('Polar Position');
      ylabel('Signal''s Intensity');

      axes_position = [0 0.05 0.5 0.95];

      ylimsize = [-1.1 1.1] * 30/2;
      xlimsize = [-1.1 1.1] * 50/2;
  end

  %hax = axes('position', axes_position);
  hax = axes;
  set(hax, ...
      'FontSize',18, ...
      'NextPlot','add', ...
      'OuterPosition', axes_position, ...
      'Visible','off', ...
      'Tag', 'axes',...
      'Parent',hfig);

  if(strncmp(type,'disp',4)) 
    himg = imshow(medfilt2(load_data(img_path,user_interf.current_frame),[5 5]));
  else
    himg = imshow(load_data(img_path,user_interf.current_frame));
  end 
  set(himg, 'Parent', hax, ...
      'Tag', 'img', ...
      'EraseMode', 'none');
  set(hax,'Ylim',ylimsize);
  set(hax,'Xlim',xlimsize);

  set(hax,'DataAspectRatio',[1 1 1]);

  hwindow.hfig = hfig;
  hwindow.hax = hax;
  hwindow.himg = himg;
  hwindow.hui = hui;
  hwindow.hgroup = hgroup;

  if (length(new_hwindows)>1)
    new_hwindows(1:end-1) = hwindows(:);
  end
  new_hwindows(end) = hwindow;

  if(strcmpi(type,'evol'))
    set(hfig,'Visible','off');
  end
  %hwindows = new_hwindows;

  return;
end

%***********************************************************
function deletecp(i)
% Delete a control point.
%
% USAGE:  deletecp(i)
%
% The control point with index i (1 to n) is deleted.
% In the case of a closed curve, a call to this function
% with i = 1 deletes both the first and last control points
% (which coincide), but a call with i = n deletes only
% the last control point and the last curve segment, thus
% converting the curve to open.

 x = mov_data.splines(user_interf.selected,user_interf.current_frame).x;
 y = mov_data.splines(user_interf.selected,user_interf.current_frame).y;
 n = hwindows(user_interf.index_main).hsplines(user_interf.selected).n;
 t = hwindows(user_interf.index_main).hsplines(user_interf.selected).t;
 hcp = hwindows(user_interf.index_main).hsplines(user_interf.selected).hcp;

% Shift arrays down.

j = i+1:n;               % indices of control points
x(j-1) = x(j); 
y(j-1) = y(j);

if i <= n-1
   delete(hcp(i));
end
j = i+1:n-1;
hcp(j-1) = hcp(j);

% Adjust tags in hcp0 and harrow.

for j = i:n-2
   set(hcp(j),'Tag',num2str(j))
end

x(n) = [];
y(n) = [];
t(n) = [];

if i <= n-1
   hcp(n-1) = [];
end

% Update n.

n = n - 1;
if i == 1
   x(n) = x(1);
   y(n) = y(1);
end

 mov_data.splines(user_interf.selected,user_interf.current_frame).x = x;
 mov_data.splines(user_interf.selected,user_interf.current_frame).y = y;
 hwindows(user_interf.index_main).hsplines(user_interf.selected).n = n;
 hwindows(user_interf.index_main).hsplines(user_interf.selected).t = t;
 hwindows(user_interf.index_main).hsplines(user_interf.selected).hcp = hcp;

return;
end  % deletecp

%***********************************************************
function dragcp
% Initiate a control point drag operation.
%
% USAGE:  dragcp;
%
% Create a dotted line (initially not visible) connecting 
% the control point with handle hc to its (one or two) 
% neighbors.  The WindowButtonMotion and WindowButtonUp 
% callbacks will update the line object and curve segments.

 x = mov_data.splines(user_interf.selected,user_interf.current_frame).x;
 y = mov_data.splines(user_interf.selected,user_interf.current_frame).y;
 n = hwindows(user_interf.index_main).hsplines(user_interf.selected).n;
 hcp = hwindows(user_interf.index_main).hsplines(user_interf.selected).hcp;
 hc = hwindows(user_interf.index_main).hc;

set(hc,'Visible','off');
set(hc,'MarkerEdgeColor',user_interf.sl_color)
set(hc,'MarkerFaceColor',user_interf.sl_color)
user_interf.drag_cp = true;

i = str2double(get(hc,'Tag'));
mp = get(hwindows(user_interf.index_main).hfig,'CurrentPoint');

xd(2) = x(i);  yd(2) = y(i);
i1 = i-1;
if i1 < 1
   i1 = n-1;
end
if i1 > 0
   xd(1) = get(hcp(i1),'XData');
   yd(1) = get(hcp(i1),'YData');
else
   xd(1) = x(i);  yd(1) = y(i);
end
i2 = i+1;
if i2 > n-1
   i2 = 1;
end
if i2 <= n-1
   xd(3) = get(hcp(i2),'XData');
   yd(3) = get(hcp(i2),'YData');
else
   xd(3) = x(i);  yd(3) = y(i);
end

% The line object consists of three points:  the control
% point and its neighbors, with the control point duplicated
% if it is an endpoint of an open curve.

hwindows(user_interf.index_main).hsplines(user_interf.selected).hdragcp = line(xd,yd,'Color',user_interf.sl_color, ...
               'LineStyle',':','Marker','o', ...
               'MarkerSize',user_interf.ms2,'Tag',num2str(i), ...
               'Visible','on');

return
end  % dragcp


%***********************************************************
function insertcp(i,p)
% Insert a control point on a curve segment or polygon edge.
%
% USAGE:  insertcp(i,p)
%
% A new control point with axes coordinates p and index i 
% (1 to n+1) is inserted (2 <= i <= n) or appended (i = 1 or
% i = n+1).  For a closed curve i is in the range 1 to n.

% Increase array lengths and shift arrays up.

 x = mov_data.splines(user_interf.selected,user_interf.current_frame).x;
 y = mov_data.splines(user_interf.selected,user_interf.current_frame).y;
 n = hwindows(user_interf.index_main).hsplines(user_interf.selected).n;
 t = hwindows(user_interf.index_main).hsplines(user_interf.selected).t;
 hcp = hwindows(user_interf.index_main).hsplines(user_interf.selected).hcp;
 hcs = hwindows(user_interf.index_main).hsplines(user_interf.selected).hcs;
 cp_color = mov_data.splines(user_interf.selected,user_interf.current_frame).cp_color;

j = i:n;                   % indices of control points
x(j+1) = x(j); 
y(j+1) = y(j);
t(j+1) = t(j);

j = i:n-1;
hcp(j+1) = hcp(j);

% Adjust tags in hcp0 and harrow.

for j = i+1:n-1+1
   set(hcp(j),'Tag',num2str(j))
end

% Increment n, n-1.

n = n + 1;

% Insert p as the new control point i.

x(i) = p(1);
y(i) = p(2);

% Select reasonable parameter values associated with the new
% control point.

if i == 1
   t(i) = t(i+1);
elseif i < n                   % 1 < i < n
   t(i) = (t(i-1)+t(i+1))/2;
else
   t(i) = t(i-1);
end
hcp(i) = line(x(i),y(i), ...
               'ButtonDownFcn',@CpButtonDnFcn, ...
               'HandleVisibility','callback', ...
               'HitTest','on', ...
               'LineStyle','none','Marker','o', ...
               'MarkerEdgeColor',cp_color, ...
               'MarkerFaceColor',cp_color, ...
               'MarkerSize',user_interf.ms2, ...
               'SelectionHighlight','off', ...
               'Tag',num2str(i));

if  strcmp(user_interf.operation,'insert_cp')
   set(hcp(i),'HitTest','off')
end
		
if strcmp(user_interf.operation,'delete_cp')  || ...
  strcmp(user_interf.operation,'move_cp')
   set(hcs,'HitTest','off')
end
	
 mov_data.splines(user_interf.selected,user_interf.current_frame).x = x;
 mov_data.splines(user_interf.selected,user_interf.current_frame).y = y;
 hwindows(user_interf.index_main).hsplines(user_interf.selected).n = n;
 hwindows(user_interf.index_main).hsplines(user_interf.selected).t = t;
 hwindows(user_interf.index_main).hsplines(user_interf.selected).hcp = hcp;
 hwindows(user_interf.index_main).hsplines(user_interf.selected).hcs = hcs;	

return;
end  % insertcp

%***********************************************************
function tsknots
% Adapt displayed knots to actual spline points

  x = mov_data.splines(user_interf.selected,user_interf.current_frame).x;
  y = mov_data.splines(user_interf.selected,user_interf.current_frame).y;
  cp_color = mov_data.splines(user_interf.selected,user_interf.current_frame).cp_color;
  if(length(cp_color)==0)
     cp_color = mov_data.splines(user_interf.selected,user_interf.previous_frame).cp_color;
     mov_data.splines(user_interf.selected,user_interf.current_frame).cp_color = cp_color;
  end
  hcp = hwindows(user_interf.index_main).hsplines(user_interf.selected).hcp;

  old_m = length(hcp);
  m = length(x);

  if(old_m>m)
    hcp(m+1:old_m) = [];
  elseif(old_m<m)

    for i=old_m+1:m
       hcp(i) = line(x(i),y(i), ...
                      'ButtonDownFcn',@CpButtonDnFcn, ...
                      'HandleVisibility','callback', ...
                      'HitTest','off', ...
                      'Visible','off', ...
                      'LineStyle','none','Marker','o', ...
                      'MarkerEdgeColor',cp_color, ...
                      'MarkerFaceColor',cp_color, ...
                      'MarkerSize',user_interf.ms2, ...
                      'Parent',hwindows(user_interf.index_main).hax, ...
                      'SelectionHighlight','off', ...
                      'Tag',num2str(i));
    end
  end

  for i=1:m
    set(hcp(i),'XData',x(i), 'YData',y(i));
  end

  hwindows(user_interf.index_main).hsplines(user_interf.selected).n = m;
  hwindows(user_interf.index_main).hsplines(user_interf.selected).hcp = hcp;

  return;
end

%***********************************************************
function tscurve
% Create and display a tension spline curve
%
% This is a nested function that provides an interface 
% between TSPACKGUI (and its callbacks) and the tension 
% spline package TSPACK. 

 x = mov_data.splines(user_interf.selected,user_interf.current_frame).x;
 y = mov_data.splines(user_interf.selected,user_interf.current_frame).y;
 if(x(end)~=x(1) | y(end)~=y(1))
   x = [x x(1)];
   y = [y y(1)];
 end

 [ t,ier] = arcl2d(x,y);
   if  ier > 0
      error('%s\n','A pair of adjacent control points coincide.')
   end

               if ier > 0
         error('Invalid knots returned by ARCLxD.')
   end

if ier == -4 
   error('A pair of adjacent control points coincide.')
end
if ier == -5
   error('Invalid bounds constraints.')
end
if ier < 0
   error('Error flag %0.0f returned by function TSPxx.\n', ier)
end

% Update the curve segment line objects.

  spline = mycsape(t,[x;y],'periodic');
  [points,distance] = fnplt(spline);
  set(hwindows(user_interf.index_main).hsplines(user_interf.selected).hcs,'XData',points(1,:), 'YData',points(2,:))

  if(length(mov_data.splines(user_interf.selected,user_interf.current_frame).cs_color)==0)
     mov_data.splines(user_interf.selected,user_interf.current_frame).cs_color = mov_data.splines(user_interf.selected,user_interf.previous_frame).cs_color;
  end
  
  hwindows(user_interf.index_main).hsplines(user_interf.selected).t = t;
  %hsplines(user_interf.selected).dist = distance;
  npixels = 5;

  for i=1:length(hwindows)
    if(i~=user_interf.index_main)
      if(i==user_interf.index_ref)
        if(~isfield(mov_data.ref_spline,'spline'))
          [t, ierr] = arcl2d(mov_data.ref_spline.x([1:end 1]), mov_data.ref_spline.y([1:end 1]));
          mov_data.ref_spline.spline = mycsape(t, [mov_data.ref_spline.x([1:end 1]); mov_data.ref_spline.y([1:end 1])]);
        end
        if(user_interf.selected==1)
          [points,distance] = fnplt(mov_data.ref_spline.spline);
          set(hwindows(i).hsplines(user_interf.selected).hcs,'XData',points(1,:), 'YData',points(2,:))
  %        set(hwindows(i).hsplines(user_interf.selected).hcs,'XData',points(1,:), 'YData',points(2,:))
        else
          [points,polar_warp] = warp_spline(mov_data.splines, mov_data.ref_spline, user_interf.selected, user_interf.current_frame, mymovie);
          set(hwindows(i).hsplines(user_interf.selected).hcs,'XData',points(:,1), 'YData',points(:,2));

          npoints = 200;
          polar_points = [0:2*pi/npoints:2*pi];
          polar_points = polar_points(1:end-1)';
          hax = findobj('Tag','profile_plot');
          hlines = get(hax, 'Children');
          hlines = hlines(end:-1:1);

          surf_profiles = zeros(length(polar_warp),length(hwindows));
          surf_colors = zeros(3,length(hwindows));

          skip = 0;
          for j=1:length(hwindows)
            if(j==user_interf.index_main || j==user_interf.index_ref)
              skip = skip + 1;
              continue;
            end

            img_path = hwindows(j).img_path;
            index = user_interf.current_frame;

            gsigma = 0.5;
            filter = fspecial('gaussian', npixels, gsigma);

            %max_dist = distance(end);
            %profile = get_spline_profile(load_data(img_path,index),spline,filter,max_dist*points);
            [profile] = get_spline_profile(load_data(img_path,index),[x' y'],filter,polar_points, mymovie, index);
            if(isfield(mymovie,'isinverted') && mymovie.isinverted)
              midindx = find(polar_points>pi,1);
            else
              profile = profile([midindx:end 1:midindx-1]);
            end
            [warp_profile] = get_spline_profile(load_data(img_path,index),[x' y'],filter,polar_warp', mymovie, index);
            
            if(length(hlines)>=j-skip)
              set(hlines(j-skip),'XData',polar_points,'YData',profile);
            else
            %profile_color = find_color(mymovie,img_path);
            profile_color = [1 0 0];

              hlines(end+1) = line(polar_points,profile, ...
                  'Color',profile_color, ...
                  'HandleVisibility', 'callback', ...
                  'HitTest','off', ...
                  'LineStyle','-', 'Marker','none', ...
                  'Parent',hax, ...
                  'SelectionHighlight','off');

              [h1,h2,hlegends,legends] = legend(hax);
              legends = cat(1,legends,get(hwindows(j).hfig,'Name'));
              hlegends = cat(1,hlegends,hlines(end));
              legend(hax, hlegends, legends, 'Location', 'NorthWest');
              legend(hax,'boxoff');
            end

            surf_colors(:,j) = get(hlines(j-skip), 'Color');
            surf_profiles(:,j) = warp_profile;

           % color_code = get(hlines(end), 'Color');
           % hsurfs = findobj('Tag','alphasurf');

          end

          intensities = max(surf_profiles,[],2);
          intensities(intensities==0) = 1e-5;
          ratios = surf_profiles ./ intensities(:,ones(1,length(hwindows)));
          cdata = zeros(length(intensities),2, 3);
          cdata(:,1,:) = ratios * surf_colors';
          cdata(:,2,:) = cdata(:,1,:);
          hsurfs = findobj('Tag',['surf' num2str(user_interf.selected)]);

          if(ishandle(hsurfs))
           set(hsurfs, 'XData', [points(:,1) points(:,1)], ...
            'YData', [points(:,2) points(:,2)], ...
            'ZData', zeros(length(points),2), ...
            'AlphaData', [intensities intensities], ...
            'CData', cdata);
          else
           hsurfs = surface('XData', [points(:,1) points(:,1)], ...
            'YData', [points(:,2) points(:,2)], ...
            'ZData', zeros(length(points),2), ...
            'Parent', hwindows(i).hax, ...
            'AlphaData', [intensities intensities], ...
            'EdgeAlpha','interp', ...
            'FaceColor', 'none', ...
            'EdgeColor', 'interp', ...
            'CData', cdata, ...
            'CDataMapping', 'direct', ...
            'Tag', ['surf' num2str(user_interf.selected)], ...
            'LineWidth', 4);
          end
        end
      else
        if(length(hwindows(i).hsplines)>=user_interf.selected)
          delete(hwindows(i).hsplines(user_interf.selected).hcs(ishandle(hwindows(i).hsplines(user_interf.selected).hcs)));
        end
        hwindows(i).hsplines(user_interf.selected).hcs = copyobj(hwindows(user_interf.index_main).hsplines(user_interf.selected).hcs,hwindows(i).hax);

        if(user_interf.selected~=1)
          [inx, iny] = get_spline_avg([x' y'], npixels, mymovie, user_interf.current_frame);

        %cdata = ones([size(surfx) 3]);
        %cdata(:,1,:) = repmat(mov_data.splines(user_interf.selected,user_interf.current_frame).cs_color,size(surfx,1),1);
        %cdata(:,2,:) = cdata(:,1,:);

        %hsurf = surface('XData', surfx, ...
        %    'YData', surfy, ...
        %    'ZData', ones(size(surfx)), ...
        %    'Parent', hwindows(i).hax, ...
        %    'CData', cdata);

          %line(surfx(:,1),surfy(:,1),'Parent',hwindows(i).hax,'Color',mov_data.splines(user_interf.selected,user_interf.current_frame).cs_color);
        htmp = line(inx,iny, ...
            'Parent',hwindows(i).hax, ...
            'Color',mov_data.splines(user_interf.selected,user_interf.current_frame).cs_color, ...
              'HandleVisibility', 'callback', ...
              'HitTest','off', ...
              'LineStyle','-', 'Marker','none', ...
              'SelectionHighlight','off');

        hwindows(i).hsplines(user_interf.selected).hcs = [htmp hwindows(i).hsplines(user_interf.selected).hcs];
        end
      end
    end
  end

  mov_data.splines(user_interf.selected,user_interf.current_frame).ppform = spline;
  mov_data.splines(user_interf.selected,user_interf.current_frame).dist = distance;

  drawnow;

  for i=1:length(hwindows)
    refresh(hwindows(i).hfig);
  end

return;
end  % tscurve

%***********************************************************
function ulpush(opcode,uindex)
% Add an entry to the undo list.
%
% USAGE:  ulpush(opcode,index)
%
% The one-character code, opcode, specifies the user_interf.operation or 
% change of curve type being performed (and to be reversed 
% if the entry is popped off the stack by Function UndoFcn),
% and index is the index of the relevant control point 
% (opcode = 'a', 'd', 'i', 'm', or 'w'), tension factor 
% (opcode = 'b' or 's'), or old curve type (opcode = 'c', 
% 'e', 'n', 't', or 'x'), where the curve type index (1:3) 
% is position in the Curve Type menu.  The following codes 
% specify user_interf.operations.
%
%    opcode = 'm':  move control point.
%    opcode = 'i':  insert control point.
%    opcode = 'd':  delete control point.
%    opcode = 'a':  alter derivative.
%    opcode = 's':  alter tension factor.
%    opcode = 'w':  alter smoothing weight.
%    opcode = 'b':  alter bounds.
%
% The following codes specify changes of curve type.
%
%    opcode = 'n/N':  ndof = 1,2,3.
%    opcode = 'c/C':  deriv = 'c2', 'c1', 'user'.
%    opcode = 't':    tension = 'shape', 'bounds', 'user'.
%    opcode = 'e':    endcond = 'periodic', 'auto', 'user'.
%    opcode = 'x/X':  interpolate, approximate.

 uptr = mov_data.splines(user_interf.selected,user_interf.current_frame).uptr;
 ulcnt = mov_data.splines(user_interf.selected,user_interf.current_frame).ulcnt;
 ulist = mov_data.splines(user_interf.selected,user_interf.current_frame).ulist;

uptr = uptr + 1;
if uptr > user_interf.upmax, uptr = 1; end
i = uindex;
data = [];

switch opcode
case {'d','m'}
   data = [mov_data.splines(user_interf.selected,user_interf.current_frame).x(i) mov_data.splines(user_interf.selected,user_interf.current_frame).y(i)];
end

ulist(uptr).opcode = opcode;
ulist(uptr).index = uindex;
ulist(uptr).data = data;
if ulcnt < user_interf.upmax
   ulcnt = ulcnt + 1;
end

 mov_data.splines(user_interf.selected,user_interf.current_frame).uptr = uptr;
 mov_data.splines(user_interf.selected,user_interf.current_frame).ulcnt = ulcnt;
 mov_data.splines(user_interf.selected,user_interf.current_frame).ulist = ulist;

return;
end  % ulpush

%***********************************************************
function AdaptSpline(hobj,event_data)

   hc = hwindows(user_interf.index_main).hc;
   if ishandle(hc)
      if strcmp(get(hc,'Marker'),'none')
         set(hc,'Color',mov_data.splines(user_interf.selected,user_interf.current_frame).cs_color)
      else
         set(hc,'MarkerEdgeColor',mov_data.splines(user_interf.selected,user_interf.current_frame).cp_color)
         set(hc,'MarkerFaceColor',mov_data.splines(user_interf.selected,user_interf.current_frame).cp_color)
      end
   end

  htmp = findobj('Tag','uiadapt');
  hwindows(user_interf.index_main).hc = [];
  if (get(htmp,'Value'))
    set(hwindows(user_interf.index_main).hui,'Enable','off');
    set(findobj('Tag','uiadapt'),'Enable','on');
    set(findobj('Tag','uimenufile'),'Enable','off');
    for i=1:mov_data.nsplines
      set(hwindows(user_interf.index_main).hsplines(i).hcs,'Hittest','on');
    end
    if(mov_data.nsplines==1)
      user_interf.operation = '';
      user_interf.selected = 1;
      set(findobj('Tag','uimenuop'),'Enable','on');
      set(hwindows(user_interf.index_main).hsplines(1).hcp,'Hittest','off','Visible','on');
    else
      user_interf.operation = 'select_sp';
      user_interf.selected = 0;
    end
  else
    set(hwindows(user_interf.index_main).hui,'Enable','on');
    set(findobj('Tag','uimenufile'),'Enable','on');
    set(findobj('Tag','uimenuop'),'Enable','off');
    for i=1:mov_data.nsplines
      set(hwindows(user_interf.index_main).hsplines(i).hcs,'Hittest','off');
      set(hwindows(user_interf.index_main).hsplines(i).hcp,'Hittest','off','Visible','off');
    end

    user_interf.operation = '';
    user_interf.selected = 0;
  end

  return;
end

%***********************************************************
function GetProfiles(hobj,event_data)

  if (mov_data.nsplines < 1)
    errordlg({'No spline can be selected.'}, ...
           'Invalid Operation')
    return;
  end

  hon = findobj('Enable','on');
  set(hon,'Enable','off');
  waitbar(0, hwait, 'Computing the splines profiles...');
  set(hwait,'Visible','on');

  npoints = 200;
  npixels = 5;
  gsigma = 0.5;
  filter = fspecial('gaussian', npixels, gsigma);

  nlayers = length(mymovie);
  subsize = ceil(sqrt(nlayers));
  init_frame = user_interf.current_frame;

  nprofiles = mov_data.nsplines*mov_data.nframes;

  for i=2:mov_data.nsplines

    points = [0:1/npoints:1];
    values = zeros(mov_data.nframes,length(points),nlayers);
    for j=1:mov_data.nframes
      user_interf.previous_frame = user_interf.current_frame;
      user_interf.current_frame = j;
      DisplayFrame;
      max_dist = mov_data.splines(i,j).dist(end);
    	
      frame_values = get_spline_profile(load_data(mymovie.maximas,j),mov_data.splines(i,j).ppform,filter,points*max_dist);
      mov_data.splines(i,j).profile = frame_values;
      values(j,:,:) = frame_values;

      waitbar(double((i-1)*mov_data.nframes + j)/nprofiles,hwait);
    end
    htmp(2*i-1) = figure;
    set(htmp(2*i-1),'Name',mov_data.splines(i,init_frame).name, ...
                'Visible', 'off','NumberTitle','off');
    for k=1:nlayers
      subplot(subsize,subsize,k);
      if(nlayers>1)
      switch k
        case 1
          subtitle = 'red channel';
        case 2
          subtitle = 'green channel';
        case 3
          subtitle = 'blue channel';
        case 4
          subtitle = 'alpha channel';
        default
          subtitle = ['channel n ' num2str(k)];
      end
      else
        subtitle = 'grayscale';
      end
      if(mov_data.nframes>1)
        %surf([0:1/npoints:1],[1:nframes],values(:,:,1)./maxint);
        surf([0:1/npoints:1],[1:mov_data.nframes],values(:,:,1));
      else
        %plot([0:1/npoints:1],(values./maxint),'Color',splines(i,user_interf.current_frame).cs_color);
        plot([0:1/npoints:1],values,'Color',mov_data.splines(i,user_interf.current_frame).cs_color);
      end
      title(subtitle);
    end

    htmp(2*i) = figure;
    set(htmp(2*i),'Name',mov_data.splines(i,init_frame).name, ...
                'Visible', 'off','NumberTitle','off');
    for k=1:nlayers
      subplot(subsize,subsize,k);
      if(nlayers>1)
      switch k
        case 1
          subtitle = 'red channel';
        case 2
          subtitle = 'green channel';
        case 3
          subtitle = 'blue channel';
        case 4
          subtitle = 'alpha channel';
        default
          subtitle = ['channel n ' num2str(k)];
      end
      else
        subtitle = 'grayscale';
      end
      if(mov_data.nframes>1)
        %imagesc([0:1/npoints:1],[1:nframes],values(:,:,1)./maxint);
        imagesc([0:1/npoints:1],[1:mov_data.nframes],values(:,:,1));
      end
      title(subtitle);
    end

  end

  %close(hwait);
  set(hwait,'Visible','off');
  set(hon,'Enable','on');
  set(htmp,'Visible','on');

  user_interf.previous_frame = user_interf.current_frame;
  user_interf.current_frame = init_frame;
  DisplayFrame;

  return;
end

%***********************************************************
function DetectElement(hobj,event_data)

  tag = get(hobj, 'Tag');
  val = get(hobj, 'Value');
  switch(tag(3:end))
    case 'cortex'
      if(val)
        if(mov_data.nsplines>1)
          set(hwindows(user_interf.index_main).hsplines(2).hcs,'Visible','on');
          DisplayFrame;
        else
          create('cortex');
        end
      else
        set(hwindows(user_interf.index_main).hsplines(2).hcs,'Visible','off');
      end
  end

  return;
end

%***********************************************************
function InvertAP(hobj,event_data)
  mymovie.isinverted = get(hobj, 'Value');

  DisplayFrame;

  return;
end
%***********************************************************
function PlayMovie(hobj,event_data)

  if(get(hobj, 'Value'))
    if (user_interf.current_frame >= mov_data.nframes)
      user_interf.previous_frame = user_interf.current_frame;
      user_interf.current_frame = 0;
    end
    start(htimer);
  else
    stop(htimer);
  end

  return;
end

%***********************************************************
function Redetect(hobj,event_data)

  hdetect = findobj('Tag','uidetectlvl');
  hfrac = findobj('Tag','uifrac');
  detect = get(hdetect,'Value');
  frac = get(hfrac,'Value');

  hon = findobj('Enable','on');
  %set(hon,'Enable','off');
  waitbar(0, hwait, 'Redetecting the splines...');
  set(hwait,'Visible','on');

  nsplines = size(mov_data.splines,1);

  for i=1:nsplines
    user_interf.selected = i;
    if(i==1)
      ppmask = [];
    else
      ppmask = mov_data.splines(1,user_interf.current_frame).ppform;
    end

    img = load_data(mymovie.edges,user_interf.current_frame);
    [x,y] = detect_border(img,ppmask,detect,false);
    if(length(x)~=0)
        [x,y,t] = optimize_spline(x,y,img,frac,0.01);

        mov_data.splines(user_interf.selected,user_interf.current_frame).x = x;
        mov_data.splines(user_interf.selected,user_interf.current_frame).y = y;
    end

    tsknots;
    tscurve;
      
    waitbar(i/nsplines,hwait);
  end

  user_interf.old_frac = frac;
  user_interf.old_lvl = detect;

  user_interf.selected = 0;
  set(hwait,'Visible','off');
  %set(hon,'Enable','on');

  return;
end

%***********************************************************
function SelectFrame(hobj,event_data)

  user_interf.previous_frame = user_interf.current_frame;
  user_interf.current_frame = round( get(hobj, 'Value' ) );
  DisplayFrame;

  return;
end

%***********************************************************
function NextFrame(hobj,event_data)

  if (user_interf.current_frame >= mov_data.nframes)
    hplay = findobj('Tag','uiplay');
    set(hplay, 'Value', 0);
    PlayMovie(hplay,event_data);
  else
    user_interf.previous_frame = user_interf.current_frame;
    user_interf.current_frame = user_interf.current_frame + 1;
    set(findobj('Tag','uiframes'), 'Value', user_interf.current_frame);
    DisplayFrame;
  end

  return;
end

%***********************************************************
function DisplayFrame

  for i=1:length(hwindows)
    if(strncmp(get(hwindows(i).hfig,'Tag'),'disp',4)) 
      set(hwindows(i).himg,'CData', medfilt2(load_data(hwindows(i).img_path,user_interf.current_frame),[5 5]));
    else
      set(hwindows(i).himg,'CData', load_data(hwindows(i).img_path,user_interf.current_frame));
    end
  end

  set(findobj('Tag','uiframestext'), 'String', sprintf('%i / %i', user_interf.current_frame, mov_data.nframes));
  
  frac = get(findobj('Tag','uifrac'),'Value');

  for i=1:size(mov_data.splines,1)
    user_interf.selected = i;
    if(i==1)
      ppmask = [];
    else
      ppmask = mov_data.splines(1,user_interf.current_frame).ppform;
    end

    if (length(mov_data.splines(i,user_interf.current_frame).x)==0)
      %if ~(get(findobj('Tag','uicopy'),'Value'))

        %img = load_data(mymovie.edges,user_interf.current_frame);
        %[x,y] = detect_border(img,ppmask,get(findobj('Tag','uidetectlvl'),'Value'),false);
        %if(length(x)~=0)
        %    [x,y,t] = optimize_spline(x,y,img,frac,0.01);
        if(i==1)
          mov_data.splines(i,user_interf.current_frame).x = mymovie.egg(user_interf.current_frame).path(:,1)';
          mov_data.splines(i,user_interf.current_frame).y = mymovie.egg(user_interf.current_frame).path(:,2)';
        else
          mov_data.splines(i,user_interf.current_frame).x = mymovie.cortex(user_interf.current_frame).path(:,1)';
          mov_data.splines(i,user_interf.current_frame).y = mymovie.cortex(user_interf.current_frame).path(:,2)';
        end
        %else
        %  if(length(mov_data.splines(i,user_interf.current_frame).x) == 0)
        %    mov_data.splines(i,user_interf.current_frame) = mov_data.splines(i,user_interf.previous_frame);
        %  end
        %end
        tsknots;
        tscurve;
      %else
      %  mov_data.splines(i,user_interf.current_frame) = mov_data.splines(i,user_interf.previous_frame);
      %end
    else
      tsknots;
      tscurve;
    end
  end
  
%  
% if(get(findobj('Tag','uicortex'),'Value'))
%  if(length(mov_data.splines(2,user_interf.current_frame).x) == 0)
%    if(get(findobj('Tag','uirecompute'),'Value'))
%
%
%        user_interf.selected = 2;
%
%        %img = get(himg,'CData');
%        %if(size(img,3)>1)
%        %  img = rgb2gray(img);
%        %end
%        img = load_data(mymovie.edges,user_interf.current_frame);
%        egg_mask = create_mask(mov_data.splines(1,user_interf.current_frame).ppform,img);
%        img = roifill(img,~egg_mask);
%        [x,y] = detect_border(img,false,get(findobj('Tag','uidetectlvl'),'Value'),false);
%        if(length(x)~=0)
%        [x,y,t] = optimize_spline(x,y,mov_data.min_dist/2);
%
%      mov_data.splines(user_interf.selected,user_interf.current_frame).x = x';
%      mov_data.splines(user_interf.selected,user_interf.current_frame).y = y';
%      
%      else
%                mov_data.splines(user_interf.selected,user_interf.current_frame).x = mov_data.splines(user_interf.selected,user_interf.previous_frame).x;
%      mov_data.splines(user_interf.selected,user_interf.current_frame).y = mov_data.splines(user_interf.selected,user_interf.previous_frame).x;
%      end
%
%        tsknots;
%        tscurve;
%
%    else
%      mov_data.splines(2,user_interf.current_frame) = mov_data.splines(2,user_interf.previous_frame);
%    end
%  else
%
%    user_interf.selected = 2;
%	tsknots;
%    tscurve;
%
% end
%end
  user_interf.selected = 0;

  for i=1:length(hwindows)
    refresh(hwindows(i).hfig);
  end
    
  return;
end

%***********************************************************
function SelectFcn(hobj,event_data)

    if(mov_data.nsplines==1)
      return;
    end

    if(user_interf.selected~=0)
    set(hwindows(user_interf.index_main).hsplines(user_interf.selected).hcp,'Hittest','off','Visible','off');
    end
    for i=1:mov_data.nsplines
      set(hwindows(user_interf.index_main).hsplines(i).hcs,'Hittest','on');
    end
    hc=[];

      user_interf.operation = 'select_sp';
      user_interf.selected = 0;

  return;
end


%***********************************************************
function FsaveFcn(hobj,event_data)
% Callback for File menu selection:  Save.

  [filename,path] = uiputfile({'*.mat'}, 'Save user_interf.current cell');
  if isequal(filename,0)  ||  isequal(path,0)
    return;
  end

  fname = fullfile(path,filename);
  save(fname, 'mov_data', 'mymovie');

  return;
end  % FsaveasFcn

%***********************************************************
function FopenFcn(hobj,event_data)
  % Callback for File menu selection:  Open.

  [filename,path] = uigetfile({'*.*'}, 'Load a file');
  if isequal(filename,0)  ||  isequal(path,0)
    return;
  end
  
  fname = fullfile(path,filename);

  switch fname(end-2:end)
    case 'mat'

    if(length(whos('-file', fname))==2)
      set(findobj('Tag','uiframes'), 'Value', user_interf.current_frame);
      new_data = load(fname, 'mov_data', 'mymovie');
      mov_data.splines = new_data.mov_data.splines;
      if(length(new_data.mymovie)>1 || ~ischar(new_data.mymovie))
        mymovie.data = store_data(new_data.mymovie,'');
      end

      [imgsize, mov_data.nframes] = size_data(mymovie.data); 

      set(hax,'Ylim',[0 imgsize(1)])
      set(hax,'Xlim',[0 imgsize(2)])

      [new_nsplines,tmp] = size(mov_data.splines);
      
      new_hsplines = repmat(struct('hcp',[],'hcs',[], ...
                 'n',0,'t',[]),1,new_nsplines);

      if(new_nsplines>mov_data.nsplines)
        for i=mov_data.nsplines+1:new_nsplines
          new_hsplines(i).hcs = line(mov_data.splines(i,1).x(1),mov_data.splines(i,1).y(1), ...
                    'ButtonDownFcn',@CsButtonDnFcn, ...
                    'Color',mov_data.splines(i,1).cs_color, ...
                    'HandleVisibility', 'callback', ...
                    'HitTest','off', ...
                    'LineStyle','-', 'Marker','none', ...
                    'Parent',hax, ...
                    'SelectionHighlight','off', ...
                    'Tag',num2str(i)); 
        end
      elseif(new_nsplines<mov_data.nsplines)
        for i=new_nsplines+1:mov_data.nsplines
          delete(hwindows(user_interf.index_main).hsplines(i).hcs);
          delete(hwindows(user_interf.index_main).hsplines(i).hcp);
        end
        hwindows(user_interf.index_main).hsplines(new_nsplines+1:end) = [];
      end

      new_hsplines(1:mov_data.nsplines) = hwindows(user_interf).hsplines;
      mov_data.nsplines = new_nsplines;
      hwindows(user_inter.index_main).hsplines = new_hsplines;

      if(get(findobj('Tag','uicortex'),'Value'))
          set(hwindows(user_inter.index_main).hsplines(2).hcs,'Visible','on');
      else
          set(hwindows(user_interf.index_main).hsplines(2).hcs,'Visible','off');
      end

        index = 1;
        for i=2:length(hwindows(user_interf.index_main).hgroup)
          if strcmp(hwindows(user_interf.index_main).hgroup(i).tag,'video')
            index = i;
          end
        end
      if(mov_data.nframes>1)
        set(hwindows(user_interf.index_main).hgroup(index).hui,'Visible','on');

          % setup frame Slider
           set(findobj('Tag','uiframes'), ...
          'Min',1, 'Max',mov_data.nframes, ...
          'SliderStep',[1/(mov_data.nframes-1) 5/(mov_data.nframes-1)], ...
          'Value', 1 );
          user_interf.current_frame = 1;
          user_interf.previous_frame = 1;
          %set(himg,'CData',[]);
          DisplayFrame;
      else
        set(hwindows(user_interf.index_main).hgroup(index).hui,'Visible','off');
          user_interf.current_frame = 1;
      end
    end

    for i=1:mov_data.nsplines
      user_interf.selected = i;
      tsknots;
      tscurve;
    end

    user_interf.selected = 0;

  otherwise
    [mymovie,ierr] = open_mymovie(fname);

    if(ierr>0)
      return
    end

    mymovie = rescale_colors(mymovie);

    [imgsize, mov_data.nframes] = size_data(mymovie.data); 

    %imgsize = size(mymovie(1).data);
    %min_dist = max(imgsize) / 15;
    %nframes = size(mymovie,2);
    user_interf.current_frame = 1; 
    user_interf.previous_frame = 1;
    mov_data.nsplines = 0;
    user_interf.selected = 0;

    mov_data.splines = [];
    %hsplines = [];

    user_interf.drag_cp = false;
    %hwindows(user_interf.index_main).hsplines(user_interf.selected).hc = [];
    %hdragcp = [];
    %hui = [];
    user_interf.operation = '';

  end

  return;
end

%***********************************************************
function FexitFcn(hobj,event_data)

  if(exist('hwindows','var'))
    for i=1:length(hwindows)
      try
        delete(findall(hwindows(i).hfig));
      catch ME
        %skip
      end
    end
  end
  if(exist('hwait','var'))
    try
      delete(hwait);
    catch ME
      %skip
    end
  end

  %if(exist('mymovie','var'))
  %  tmpfiles = fieldnames(mymovie);
  %  for i=1:length(tmpfiles)
  %    delete(getfield(mymovie,char(tmpfiles{i})));
  %  end
  %end

  clear all;

  return;
end

%***********************************************************
function FvoidFcn(hobj,event_data)
  return;
end

%***********************************************************
function ButtonUpFcn(hobj,event_data)
% Callback for mouse button release.

if user_interf.drag_cp

% Button release following a control point drag:  update 
% x, y, z, the curve segments, and the line objects with
% the new position.

   user_interf.drag_cp = false;
   set(hwindows(user_interf.index_main).hc,'Visible','on');
   %set(hc,'MarkerEdgeColor',splines(user_interf.selected,user_interf.current_frame).cp_color)
   %set(hc,'MarkerFaceColor',splines(user_interf.selected,user_interf.current_frame).cp_color)
   xd = get(hwindows(user_interf.index_main).hsplines(user_interf.selected).hdragcp,'XData');
   yd = get(hwindows(user_interf.index_main).hsplines(user_interf.selected).hdragcp,'YData');
   i = str2double(get(hwindows(user_interf.index_main).hsplines(user_interf.selected).hdragcp,'Tag'));
   delete(hwindows(user_interf.index_main).hsplines(user_interf.selected).hdragcp)

% Compute the distance d in window coordinates from the new
% control point (xd(2),yd(2),zd(2)) to each of its two
% neighbors (xd(1),yd(1),zd(1)) and (xd(3),yd(3),zd(3)).
% If d < dtol then move the control point away from its
% neighbor in the case of an interpolatory or smoothing 
% curve, or create a duplicate control point in the case of
% a B-spline curve.
 
   dtol = 4;
   dely = get(hwindows(user_interf.index_main).hax,'YLim');
   fp = get(hwindows(user_interf.index_main).hfig,'Position');
   dely = 2*dtol*(dely(2)-dely(1))/fp(4);
   p = [xd(2) yd(2)];
   w = xform(p);
         pn = [xd(1) yd(1)];
         wn = xform(pn);
         d = sqrt((w(1)-wn(1))^2 + (w(2)-wn(2))^2);
      if d < dtol
            yd(2) = yd(2) + dely;
      end
         pn = [xd(3) yd(3)];
         wn = xform(pn);
         d = sqrt((w(1)-wn(1))^2 + (w(2)-wn(2))^2);
      if d < dtol
            yd(2) = yd(2) + dely;
      end

 x = mov_data.splines(user_interf.selected,user_interf.current_frame).x;
 y = mov_data.splines(user_interf.selected,user_interf.current_frame).y;

   x(i) = xd(2);
   y(i) = yd(2);
   if i == 1
      n = hwindows(user_interf.index_main).hsplines(user_interf.selected).n;
      x(n) = x(1);
      y(n) = y(1);
   end
 set(hwindows(user_interf.index_main).hsplines(user_interf.selected).hcp(i),'XData',x(i), 'YData',y(i))
 mov_data.splines(user_interf.selected,user_interf.current_frame).x = x;
 mov_data.splines(user_interf.selected,user_interf.current_frame).y = y;

   tscurve;

end

return;
end  % ButtonUpFcn


%***********************************************************
function CpButtonDnFcn(hobj,event_data)
% Callback for mouse button press on a control point.

% Set s to the selection (mouse click) type:  

s = get(hwindows(user_interf.index_main).hfig,'SelectionType');
if strcmp(s,'normal')

% Left button click:
%
% Change the user_interf.currently user_interf.selected curve object to hobj:
  hc = hwindows(user_interf.index_main).hc;
   if ishandle(hc)
      if strcmp(get(hc,'Marker'),'none')
         set(hc,'Color',mov_data.splines(user_interf.selected,user_interf.current_frame).cs_color)
      else
         set(hc,'MarkerEdgeColor',mov_data.splines(user_interf.selected,user_interf.current_frame).cp_color)
         set(hc,'MarkerFaceColor',mov_data.splines(user_interf.selected,user_interf.current_frame).cp_color)
         %set(hc,'Visible','off')
      end
   end
   hc = hobj;
   hwindows(user_interf.index_main).hc = hobj;
   i = str2double(get(hc,'Tag'));
   switch user_interf.operation

   case 'delete_cp'

% Delete control point i, and update the curve.

      if hwindows(user_interf.index_main).hsplines(user_interf.selected).n <= 5 
          errordlg({'A closed curve requires at least four ', ...
            'control points.'}, ...
           'Invalid Operation')
           user_interf.operation = '';
      else 
        ulpush('d',i)
        deletecp(i)
        tscurve;
      end

   case 'move_cp'

% Allow the user to drag the control point:  create a dotted
% line (initially not visible) connecting the point to its 
% (one or two) neighbors.  The WindowButtonMotion and
% WindowButtonUp callbacks will update the line object and 
% curve segments.
 
      ulpush('m',i)
      dragcp;
   end
elseif strcmp(s,'alt')

% Right button click:

end

return;
end  % CpButtonDnFcn

%***********************************************************
function CsButtonDnFcn(hobj,event_data)
% Callback for mouse button press on a curve segment.

% Set s to the selection type:  

s = get(hwindows(user_interf.index_main).hfig,'SelectionType');
if strcmp(s,'normal')

% Left button click:
%
% Change the user_interf.currently user_interf.selected curve object to hobj.
  hc = hwindows(user_interf.index_main).hc;
   if ishandle(hc)
      if strcmp(get(hc,'Marker'),'none')
         set(hc,'Color',mov_data.splines(user_interf.selected,user_interf.current_frame).cs_color)
      else
         set(hc,'MarkerEdgeColor',mov_data.splines(user_interf.selected,user_interf.current_frame).cp_color)
         set(hc,'MarkerFaceColor',mov_data.splines(user_interf.selected,user_interf.current_frame).cp_color)
      end
   end
   hwindows(user_interf.index_main).hc = hobj;
   switch user_interf.operation
   case 'select_sp'

    user_interf.selected = str2double(get(hobj,'Tag'));

    for i=1:mov_data.nsplines
      if (i==user_interf.selected)
        set(hwindows(user_interf.index_main).hsplines(i).hcp,'Hittest','on','Visible','on');
      else
        set(hwindows(user_interf.index_main).hsplines(i).hcs,'Hittest','off');
        set(hwindows(user_interf.index_main).hsplines(i).hcp,'Hittest','off','Visible','off');
      end
    end
    user_interf.operation = '';

    set(findobj('Tag','uimenuop'),'Enable','on');
   case 'insert_cp'

% Insert a new control point at the mouse position (which 
% lies on the user_interf.selected curve segment).

% Convert the ne points on curve segment hcs0(i-1) to window
% coordinates wp (ne by 2), and find the nearest point p to
% the mouse position mp.

      vx = get(hwindows(user_interf.index_main).hsplines(user_interf.selected).hcs,'XData');
      vy = get(hwindows(user_interf.index_main).hsplines(user_interf.selected).hcs,'YData');
      wp = [vx' vy'];
      mp = vform(get(hwindows(user_interf.index_main).hfig,'CurrentPoint'));
      wp = [wp(:,1)-mp(1) wp(:,2)-mp(2)]';
      [m,j] = min(max(abs(wp)));
      p(1) = vx(j(1));
      p(2) = vy(j(1));
      i = min(find(hwindows(user_interf.index_main).hsplines(user_interf.selected).t > mov_data.splines(user_interf.selected,user_interf.current_frame).dist(j)));
      ulpush('i',i)
      insertcp(i,p)
      hwindows(user_interf.index_main).hc = hwindows(user_interf.index_main).hsplines(user_interf.selected).hcp(i);

% Initiate a drag.

      dragcp;
   end
elseif strcmp(s,'alt')

% Right button click:

end

return;
end  % CsButtonDnFcn

%***********************************************************
function DeleteCpFcn(hobj,event_data)
% Callback for Operations menu selection:  delete control point.

if hwindows(user_interf.index_main).hsplines(user_interf.selected).n <= 4
   errordlg({'A closed curve requires at least four ', ...
            'control points.'}, ...
           'Invalid Operation')
   return;
end

user_interf.operation = 'delete_cp';
set(hwindows(user_interf.index_main).hsplines(user_interf.selected).hcp,'HitTest','on')
set(hwindows(user_interf.index_main).hsplines(user_interf.selected).hcs,'HitTest','off')
set(hwindows(user_interf.index_main).hsplines(user_interf.selected).hcp,'Visible','on')

return;
end  % DeleteCpFcn


%***********************************************************
function InsertCpFcn(hobj,event_data)
% Callback for Operations menu selection:  insert control point.

user_interf.operation = 'insert_cp';
set(hwindows(user_interf.index_main).hsplines(user_interf.selected).hcp,'HitTest','off')
set(hwindows(user_interf.index_main).hsplines(user_interf.selected).hcs,'HitTest','on')
set(hwindows(user_interf.index_main).hsplines(user_interf.selected).hcp,'Visible','on')

return;
end  % InsertCpFcn


%***********************************************************
function KeyPressFcn(hobj,event_data)
% Callback for key press in figure window.

% Get the ASCII value key.

key = double(event_data.Character);

% Shift, Alt, and Ctrl keys, when held down, generate 
% callbacks at the keyboard repeat rate, and these keys
% by themselves are stored as ''.

if isempty(key), return; end   

modifier = event_data.Modifier;
if length(modifier) > 0  &&  strcmp(modifier{1},'control')
   axes(hax);  % Return control to figure 1 for menu hotkeys.
   return;
end

%switch key
%case 28

% Left Arrow keypress:  rotate camera clockwise (as viewed
% from positive z).

%case 29

% Right Arrow keypress:  rotate camera ccw.

%case 30

% Up Arrow keypress:  increase elevation.

%case 31

% Down Arrow keypress:  decrease elevation.

%case double('z')

% z key:  Zoom in.
   
%case double('Z')

% Z key:  Zoom out.

%end

return;
end  % KeyPressFcn


%***********************************************************
function MotionFcn(hobj,event_data)
% Callback for mouse motion with a button pressed in figure 1 or figure 4.

if user_interf.drag_cp

% Drag a control point (line object with handle hdragcp
% created by CpButtonDnFcn):  convert the mouse position to
% axes coordinates and update the second (and possibly the
% first or third) point of the line object. 

   mp = get(hwindows(user_interf.index_main).hfig,'CurrentPoint');
   xd = get(hwindows(user_interf.index_main).hsplines(user_interf.selected).hdragcp,'XData'); 
   yd = get(hwindows(user_interf.index_main).hsplines(user_interf.selected).hdragcp,'YData'); 
   p = vform(mp);
   if xd(1) == xd(2)  &&  yd(1) == yd(2) 
      xd(1) = p(1);
      yd(1) = p(2);
   end
   if xd(3) == xd(2)  &&  yd(3) == yd(2)
      xd(3) = p(1);
      yd(3) = p(2);
   end
   xd(2) = p(1);
   yd(2) = p(2);
   set(hwindows(user_interf.index_main).hsplines(user_interf.selected).hdragcp,'XData',xd, 'YData',yd, ...
       'Visible','on')

%elseif hslider4 == get(hfig4,'CurrentObject')

end

return;
end  % MotionFcn

%***********************************************************
function FitSpline(hobj,event_data)

  % Find correct window :
  hfig = get(hobj,'Parent');
  himg = findobj(hfig,'Tag','img');
  img = get(himg,'CData');

  set(findobj(findobj('Tag','evol'),'Tag','img'),'CData',img);

  set(findobj('Tag','uiscoretxt'),'String',' ');
  user_interf.selected = 1;
  frac = get(findobj('Tag','uifrac'),'Value');

  % Get from popupmenu item :
  %method = get(popupmenu,String){get(popupmenu,Value)}

  optimize = 'maximas';
  method = 'montecarlo';
  scoring = 'global';
  range.x = [-5:0.5:5];
  range.y = [-5:0.5:5];
  range.n  = str2num(get(findobj('Tag','uitemp'),'String'));
  converge = true;
  %[x,y,score] = fit_spline(load_data(mymovie.maximas,user_interf.current_frame),mov_data.splines(1,user_interf.current_frame).x,mov_data.splines(1,user_interf.current_frame).y,method,range,scoring,optimize,converge);
  [x,y,score] = fit_spline(img,mov_data.splines(1,user_interf.current_frame).x,mov_data.splines(1,user_interf.current_frame).y,method,range,scoring,frac,optimize,converge);
  mov_data.splines(1,user_interf.current_frame).x = x;
  mov_data.splines(1,user_interf.current_frame).y = y;

  tsknots;
  tscurve;

  set(findobj('Tag','uiscoretxt'),'String',num2str(score));
end

%***********************************************************
function Score(hobj,event_data)

 length = mov_data.splines(1,user_interf.current_frame).dist(end);
 %npoints = ceil(length/100)*100;
 npoints = 300;
 points = [0:length/(npoints-1):length];
 frac = get(findobj('Tag','uifrac'),'Value');
 
 method = 'gaussian';
 interp = 'subpixsum';
 [score, parts] = get_spline_score(mov_data.splines(1,user_interf.current_frame).ppform,load_data(mymovie.maximas,user_interf.current_frame),frac,points,method,interp);

 set(findobj('Tag','uiscoretxt'),'String',[num2str(score) ':' num2str(parts)]);

return;
end

%***********************************************************
function MoveCpFcn(hobj,event_data)
% Callback for Operations menu selection:  move control point.

user_interf.operation = 'move_cp';
set(hwindows(user_interf.index_main).hsplines(user_interf.selected).hcp,'HitTest','on')
set(hwindows(user_interf.index_main).hsplines(user_interf.selected).hcs,'HitTest','off')
set(hwindows(user_interf.index_main).hsplines(user_interf.selected).hcp,'Visible','on')

return;
end  % MoveCpFcn

%***********************************************************
function UndoFcn(hobj,event_data)
% Callback for Undo Last Operation menu selection.
%
% Refer to Function ulpush for a description of the stack
% data structure.

 uptr = mov_data.splines(user_interf.selected,user_interf.current_frame).uptr;
 ulcnt = mov_data.splines(user_interf.selected,user_interf.current_frame).ulcnt;
 ulist = mov_data.splines(user_interf.selected,user_interf.current_frame).ulist;

if ulcnt == 0
   errordlg('The previous operation cannot be undone', ...
           'Empty Stack')
   return;
end
op = ulist(uptr).opcode;
i = ulist(uptr).index;
data = ulist(uptr).data;
ulcnt = ulcnt - 1;
uptr = uptr - 1;
if uptr < 1
   uptr = user_interf.upmax;
end

switch op
case 'd'

% Re-insert the deleted control point.

   cp = [data(1) data(2)];
   insertcp(i,cp)

case 'i'

% Delete the inserted control point.

   deletecp(i)

case 'm'

% Restore the moved control point.

 x = mov_data.splines(user_interf.selected,user_interf.current_frame).x;
 y = mov_data.splines(user_interf.selected,user_interf.current_frame).y;

   x(i) = data(1);
   y(i) = data(2);
   if i == 1
      n = hwindows(user_interf.index_main).hsplines(user_interf.selected).n;
      x(n) = x(1);
      y(n) = y(1);
   end
   set(hwindows(user_interf.index_main).hsplines(user_interf.selected).hcp(i),'XData',x(i), 'YData',y(i))

 mov_data.splines(user_interf.selected,user_interf.current_frame).x = x;
 mov_data.splines(user_interf.selected,user_interf.current_frame).y = y;
end

% Compute the new curve.

tscurve;

 mov_data.splines(user_interf.selected,user_interf.current_frame).uptr = uptr;
 mov_data.splines(user_interf.selected,user_interf.current_frame).ulcnt = ulcnt;
 mov_data.splines(user_interf.selected,user_interf.current_frame).ulist = ulist;

return;
end  % UndoFcn
end  % tspackgui

%***********************************************************
%***********************************************************


%***********************************************************
function vout = vform(p)
% Convert a mouse position (window coordinates) into a vertex
% position (data space).

hax = gca;

% Compute the translation vector defining operator T1
% (first transformation) and the product xf of the 
% additional transformations.

xl = get(hax,'XLim');
yl = get(hax,'YLim');
zl = [0 1];

% Sometimes, when stretch-to-fill mode is in effect, T2(3,3) = 
% 1/(zl(2)-zl(1)) is equal to 5.0e15, causing xf to be nearly
% singular.

xf = get(hax,'x_RenderTransform');
if abs(xf(3,3)) > 5.0e15
   xf(3,3) = xf(3,3)/5.0e15;
end

% Apply the translation operator:  vt = T1*v.

vt = zeros(3,1);
vt(1) = - xl(1);
vt(2) = - yl(1);
vt(3) = - zl(1);

% Compute the depth associated with v:  3rd component of xf*vt.

depth = xf(3,1)*vt(1) + xf(3,2)*vt(2) + ...
        xf(3,3)*vt(3) + xf(3,4);

% Convert the mouse position to device coordinates relative
% to the upper left corner of the user_interf.current figure window,
% solve xf*vt = b for b = [p(1);p(2);depth;1], and remove
% the fourth component (which should be 1).

b = get(gcf,'Position');
b = [p(1); b(4)+0.5-p(2); depth; 1];
vt = xf\b;
vt(4) = [];

% Compute vout = T1^(-1)*vt.

vout = vt + [xl(1); yl(2); zl(1)];

limits = [xl; yl; zl];
vout(vout<limits(:,1)) = limits(vout<limits(:,1),1);
vout(vout>limits(:,2)) = limits(vout>limits(:,2),2);

return;

end

%***********************************************************
function [p] = xform(v)
% Transform vertices from data space to window coordinates.
%
% USAGE:  [p] = xform(v);
%
% Given an N by 3 array V containing a set of vertices (axes
% coordinates), this function returns an N by 2 array P 
% containing the corresponding window coordinates, and,
% optionally, the column vector of depths.  The window
% coordinate space is [0,w] X [0,h] for a figure window with
% width w and height h in pixels with (0,0) at the lower
% left corner.  If the vertices are points user_interf.selected by the
% mouse (CurrentPoint property of the axes), the window
% coordinates will be pixel centers of the form i-0.5 for an
% integer i in the range 1 to w or 1 to h.
%
% Note that the transformation uses properties of the 
% user_interf.current axes and the user_interf.current figure, and the Units 
% property of the figure is assumed to be pixels (the 
% default).

hax = gca;
xl = get(hax,'XLim');
yl = get(hax,'YLim');
zl = [0 1];

% Sometimes, when stretch-to-fill mode is in effect, T2(3,3) = 
% 1/(zl(2)-zl(1)) is equal to 5.0e15, causing xf to be nearly
% singular.

xf = get(hax,'x_RenderTransform');
if abs(xf(3,3)) > 5.0e15
   xf(3,3) = xf(3,3)/5.0e15;
end

% Apply the translation operator:  vt = T1*v.

vt = zeros(size(v));
vt(:,1) = v(:,1) - xl(1);
vt(:,2) = v(:,2) - yl(1);
vt(:,3) = - zl(1);

% Apply xf to the rows of vt (converted to homogeneous
% coordinates) to get (p(:,1);p(:,2);depth(:);w(:)).
% Perspective projection results in w ~= 1, requiring
% normalization.

w = xf(4,1)*vt(:,1) + xf(4,2)*vt(:,2) + xf(4,3)*vt(:,3) + xf(4,4);
p(:,1) = (xf(1,1)*vt(:,1) + xf(1,2)*vt(:,2) + ...
          xf(1,3)*vt(:,3) + xf(1,4))./w;
p(:,2) = (xf(2,1)*vt(:,1) + xf(2,2)*vt(:,2) + ...
          xf(2,3)*vt(:,3) + xf(2,4))./w;

% The window coordinate system of the transformed vertices 
% has origin at the upper left.  The y-components must be
% subtracted from the window height.

fp = get(gcf,'Position');
p(:,2) = fp(4) - p(:,2);

return;
end  % xform


%***********************************************************
% TSPACK Functions
%***********************************************************
function [t,ier] = arcl2d(x,y)
% arcl2d:  Computes cumulative arc lengths along a planar curve
%
% USAGE:  [t,ier] = arcl2d(x,y);
%
%   Given an ordered sequence of N points (X,Y) defining a
% polygonal curve in the plane, this function computes the
% sequence T of cumulative arc lengths along the curve:
% T(1) = 0 and, for 2 <= K <= N, T(K) is the sum of
% Euclidean distances between (X(I-1),Y(I-1)) and (X(I),Y(I))
% for I = 2 to K.  A closed curve corresponds to X(1) =
% X(N) and Y(1) = Y(N), and more generally, duplicate points
% are permitted but must not be adjacent.  Thus, T contains
% a strictly increasing sequence of values which may be used
% as parameters for fitting a smooth curve to the sequence
% of points.
%
% On input:
%
%       X,Y = Vectors of length N containing the coordinates
%             of the points.
%
% On output:
%
%       T = Vector of size(X) containing cumulative arc 
%           lengths defined above.
%
%       IER = Error indicator:
%             IER = 0 if no errors were encountered.
%             IER = I if X(I) = X(I+1) and Y(I) = Y(I+1) for
%                     some I in the range 1 to N-1.
%
% Modules required by ARCL2D:  None
%
%***********************************************************

% Set ds to the vector of arc lengths, and compute t.

ds = sqrt(diff(x).^2 + diff(y).^2);
m = size(x);
if (m(1) == 1)
   t = [0 cumsum(ds)];
else
   t = [0; cumsum(ds)];
end
if (nargout > 1) 

% Test for a zero arc length.

   if (all(ds))
      ier = 0;
   else
      ier = find(~ds,1);  % ier = index of first zero
   end
end
return;

end  % arcl2d

%***********************************************************
function [mymovie,ierr] = open_mymovie(fname)

  ierr = 0;

  %[file,path,filter]= uigetfile('*.*','Select a Movie or a Picture');
  %if(isequal(file,0))
  %  mymovie = [];
  %  ierr = 1;
  %  return;
  %end

  %mymovie.filename=file;
  %fname = fullfile(path,file);

  mymovie.data = '';
  switch fname(end-2:end)
  
    case 'avi'
      m=aviread(fname);

      nframes = length(m);
      %if(frames > 10)
      %  frames = 10;
      %end
      
      mymovie.data = '';
      for n=1:nframes
        mymovie.data=store_data(m(n).cdata,mymovie.data);
      end;
 
    case 'jpg'
      mymovie.data=store_data(double(imread(fname,'jpeg')),mymovie.data);

    %case {'tif','stk','lsm'}
    %  pics = tiffread(fname);
    %  for i=1:length(pics)
    %    mymovie.data = store_data(pics(i).data,mymovie.data);
    %  end

    otherwise
      %mymovie=read_movie(fname);
      mymovie = bfopen(fname);
  end

  mymovie = preprocess_movie(mymovie);
  return;
end

%***********************************************************
function [mymovie] = compute_edges(mymovie)
  
  hwait = findobj('Tag','wait');
  waitbar(0, hwait, 'Computing Edges...');
  set(hwait,'Visible','on');
  %hwait = waitbar(0,'Computing Edges...','Name','CellCoord Info');

  if(~isfield(mymovie,'edges'))
    mymovie.edges='';
  end
  if(~isfield(mymovie,'maximas'))
    mymovie.maximas='';
  end
  if(~isfield(mymovie,'angles'))
    mymovie.angles='';
  end

  [imgsize, nframes] = size_data(mymovie.dic);
  %[e,f] = size(mymovie);
  sobel = fspecial('sobel');
  lofg = fspecial('log',[5 5],0.5);
  for i=1:nframes
    img = load_data(mymovie.dic,i);
    sob1 = imfilter(img,sobel,'symmetric');
    sob2 = imfilter(img,sobel','symmetric');
    img = sqrt((sob1.^2)+(sob2.^2));
    img = img - min(min(img));
    img = img./max(max(img));
    mymovie.edges = modify_data(mymovie.edges,img,i);
    img = atan2(sob2, sob1);
    mymovie.angles = modify_data(mymovie.angles,img,i);
    img = imfilter(img,lofg,'symmetric');
    %tmpimg = tmpimg - min(min(tmpimg));
    img = img - min(min(img));
    img = img./max(max(img));
    mymovie.maximas = modify_data(mymovie.maximas,img,i);

    waitbar(i/nframes,hwait);
  end
  clear img sob1 sob2;

  %close(hwait);
  set(hwait,'Visible','off');
end

%***********************************************************
function [x,y] = detect_border(imd, ppmask, lvl, interactive)

  %if(size(imd,3)>1)
  %  imd = rgb2gray(imd);
  %end
  %[junk threshold] = edge(imd,'sobel');
  %imd = edge(imd,'sobel', threshold*lvl);

  border_size = 20;

  if(isstruct(ppmask))
    egg_mask = create_mask(ppmask,imd);
    egg_mask = imerode(egg_mask,strel('disk',3));
    imd = imd.*egg_mask;
  end

  thresh = graythresh(imd);
  a = 0.5 - thresh;
  b = 0.5 - 2*a;

  imd = (imd>((a*lvl^2) + (b*lvl)));
  imd = increase_canevas(imd,border_size);
  se90 = strel('line', 3, 90);
  se0 = strel('line', 3, 0);
  seround = strel('diamond', 4);
  imd = imdilate(imd, [se90 se0]);
  imd = imfill(imd, 'holes');
  imd = imerode(imd, seround);
  imd = reduce_canevas(imd,border_size);
  imd = bwareaopen(imd,5000);
  if(any(any(imd==1)))
  
      bounds = bwboundaries(imd,8,'noholes');
      bounds = bounds{1};
      %if (is_convex)
      if(~isstruct(ppmask))
        indexes = convhull(bounds(:,1),bounds(:,2));
        bounds = bounds(indexes,:);
      end
      %else
        x = bounds(:,2);
        y = bounds(:,1);
      %end
      %x = [x; x(1)];
      %y = [y; y(1)];
  elseif(interactive)
     warndlg('No border was detected, please create one','CellCoord Warning','modal');
      
     mainwindow = findobj('Tag','main');
     htmp = impoly(findobj(mainwindow,'Tag','axes'),[]);
     api = iptgetapi(htmp);
     pos = api.getPosition();
     %x = [pos(:,1); pos(1,1)];
     %y = [pos(:,2); pos(1,2)];
     x = pos(:,1);
     y = pos(:,2);
     delete(htmp);
  else
     %warndlg('No border was detected, keeping the old one','CellCoord','modal');
     x = [];
     y = [];
  end

  x = x';
  y = y';

  return;
end

%***********************************************************
function resized = increase_canevas(img, new_size)

  [h,w] = size(img);
  
  resized =  zeros(h+(2*new_size), w+(2*new_size));
  resized((new_size+1):(new_size+h),(new_size+1):(new_size+w)) = img;

end

%***********************************************************
function resized = reduce_canevas(img, new_size)

  [h,w] = size(img);
  
  resized =  img((new_size+1):(h-new_size), (new_size+1):(w-new_size));

end

%***********************************************************
function [bestx,besty,bestscore] = fit_spline(img,x,y,method,range,scoring,frac,optimize,converge)
  
  %Monte Carlo Fitting:
  %x0 = [0 0];
  %lb = [-64 -64];
  %ub = [64 64];
  %[x,fval] = simulannealbnd(@dejong5fcn,x0,lb,u%

  hwait = findobj('Tag','wait');
  waitbar(0, hwait, 'Fitting Spline...');
  set(hwait,'Visible','on');
  %hwait = waitbar(0,'Fitting Spline...','Name','CellCoord Info');

  %fitmov = avifile('myfitmovie');

  if nargin<9
    converge = false;
  end
  if nargin<8
    optimize = 'maximas';
  end
  if nargin<7
    frac = 0.5;
  end
  if nargin<6
    converge = 'global';
  end

  if length(x)<2
    return;
  end

  nknots = length(x)-1;

  bestx = x(1:nknots);
  besty = y(1:nknots);
  bestscore = -Inf;

  switch method
    case 'radial'
      npts = length(range.x);

      [t,ier] = arcl2d(x,y);
      spline = csape(t,[x;y],'periodic');

      deriv = fnval(fnder(spline),t);
      orth = [deriv(2,:); -deriv(1,:)];
      norm = sqrt(sum(orth.^2,1));
      orth = orth ./ [norm;norm];
  
      newx = range.x'*orth(1,:);
      newx = newx + (ones(npts,1)*x);

      newy = range.x'*orth(2,:);
      newy = newy + (ones(npts,1)*y);
    case 'neighborhood'
      nptsx = length(range.x);
      nptsy = length(range.y);
      newx = range.x'*ones(1,nknots)+ones(nptsx,1)*x(1:nknots);
      newy = range.y'*ones(1,nknots)+ones(nptsy,1)*y(1:nknots);

      newx = repmat(newx, nptsy, 1);
      newy = repmat(newy, 1, nptsx)';
      newy = reshape(newy, nknots, [])';
    %case 'symbolic'
    %  symbspline = symbolic_csape(nknots);
    case 'montecarlo'
      x0 = [x(1:nknots);y(1:nknots)];
      %lb = x0 - range.n;
      %ub = x0 + range.n;
      lb = zeros(size(2,nknots));
      ub = [size(img,2)*ones(1,nknots);size(img,1)*ones(1,nknots)];
      options = struct('FunctionParameters', struct('img', img, 'frac', frac),... %, 'fitmov', fitmov), ...
                       'AnnealingFcn', @brownianannealing, ...
                       'PlotFcns',@spline_evol, ...
                       'MaxIter', 1000, ...
                       'InitialTemperature', range.n);
                       %'Display', 'iter', ...
      [x,bestscore] = mysimulannealbnd(@get_spline_score,x0,lb,ub,options);
      x = abs(x);
      bestx = [x(1,:) x(1,1)];
      besty = [x(2,:) x(2,1)];

      %close(hwait);
      set(hwait,'Visible','off');
      return;
  end

 imgsize = size(img);
 newx((newx < 1) | (newx > imgsize(2))) = NaN;
 newy((newy < 1) | (newy > imgsize(1))) = NaN;
  
  switch optimize
    case 'none'
    case 'maximas'
      for i=1:nknots
        indx = sub2ind(size(img),round(newy(:,i)),round(newx(:,i)));
        nanindx = isnan(indx);
        vals = img(indx(~nanindx));
        [vals, indx] = sort(vals,'descend');
        indx = [indx; find(nanindx)];
        indx = indx(1:range.n);
        newx(1:range.n,i) = newx(indx,i);
        newy(1:range.n,i) = newy(indx,i);
      end
      
      newx = newx(1:range.n,:);
      newy = newy(1:range.n,:);
  end

  [npts,nknots] = size(newx);

  switch scoring
    case 'global'

      for i=1:nknots
        bestpos = 0;
        bestscore = -Inf;
        for j=1:npts

          if(isnan(newx(j,i)) || isnan(newy(j,i)))
            continue;
          end

          bestx(i) = newx(j,i);
          besty(i) = newy(j,i);

          %[ttmp,ier] = arcl2d([bestx bestx(1)],[besty besty(1)]);
          %spline = csape(ttmp,[bestx bestx(1);besty besty(1)],'periodic');
          %[points,distance] = fnplt(spline);
          %tmpscore = get_spline_score(img,spline,[0:distance(end)/1000:distance(end)],'direct','subpixsum');
          tmpscore = get_spline_score([bestx;besty],img,frac);

          if(tmpscore>bestscore)
            bestpos = j;
            bestscore = tmpscore;
          end
          waitbar(((i-1)*npts + j)/(nknots*npts),hwait);
        end
        bestx(i) = newx(bestpos,i);
        besty(i) = newy(bestpos,i);
      end
      bestx = [bestx bestx(1)];
      besty = [besty besty(1)];
  end

  %close(hwait);
  set(hwait,'Visible','off');

  if(converge)
    success = false;
    for i=1:20
      [newx,newy,newscore] = fit_spline(img,bestx,besty,method,range,scoring,frac,optimize,false);
      if(newscore<bestscore || all(newx==bestx) && all(newy==besty))
        success = true;
        break;
      end
      bestx = newx;
      besty = newy;
      bestscore = newscore;
    end

    if(~success)
     warndlg('Sorry, no convergeance was reached','CellCoord Warning','modal');
    end
  end
end
%***********************************************************
function [final_score,parts] = get_spline_score(spline,img,frac,points,method,interp)

    if(~isstruct(spline))
      [m,n] = size(spline);

      if(n<m)
        spline=spline';
        [m,n] = size(spline);
      end

      if(m==2)
        x = abs([spline(1,:) spline(1,1)]);
        y = abs([spline(2,:) spline(2,1)]);
      else
        x = [spline(1,1:end/2) spline(1,1)];
        y = [spline(1,((end/2)+1):end) spline(1,(end/2)+1)];
      end
      
      [ttmp,ier] = arcl2d(x,y);
      spline = csape(ttmp,[x;y],'periodic');
      [pttmp,distance] = fnplt(spline);

      if(isstruct(img))
        args = img;
        args_names = fields(img);
        for i=1:length(args_names)
          eval([args_names{i} '=args.' args_names{i} ';']);
        end
      end
    end
    if ~exist('interp','var')
      interp = 'subpixsum';
    end
    if ~exist('method','var')
      method = 'gaussian';
    end
    if ~exist('points','var')
      points = [0:distance(end)/1000:distance(end)];
    else
      if(length(points)==1)
        points = [0:distance(end)/points:distance(end)];
      end
    end
    if ~exist('frac','var')
      frac = 0.5;
    end

    position = get_spline_positions(spline,points,size(img));
    center = mean(position,2);
    %rad = mean(sqrt(sum((position - center).^2,1)))
    index = ~(isnan(position(1,:)));
    position = position(:,index);
    points = points(index);

    %acc = sum((fnval(fnder(spline,2),points)).^2,1);
    %acc = (fnval(fnder(spline,2),points));
    %acc = sum(sqrt(sum(acc.^2,1)));
    area = sum(fnval(fnder(spline,-2),points([1 end])),1);
    area = diff(area);
    acc = (diff(points([1 end])).^2)/(4*pi*area);
    %acc = rad*diff(points([1 end]))/area
    %acc = acc * area / 1e9;
    %acc = diff(sqrt(acc));
    %acc = sum(acc)
    %acc = sum(exp(acc))/length(points);

    %position = fnval(spline,points);

    switch method
      case 'direct'
      case 'gaussian'
        img = imfilter(img,fspecial('gaussian'),'symmetric');
      case 'log'
        img = imfilter(img,fspecial('log'),'symmetric');
      case 'invlog'
        img = imfilter(img,(-fspecial('log')),'symmetric');
    end
    switch interp
      case 'pixsum'
        position = round(position);
        ssize = size(img);
        indx = sub2ind(size(img),position(2,:),position(1,:));
        score = img(indx);
        score = sum(score) / length(points);
      case 'subpixsum'
        img = [img(1,:); img; img(end,:)];
        img = [img(:,1) img img(:,end)];
        ssize = size(img);
        position = position + 1;
        positionc = ceil(position);
        positionf = floor(position);
        distc = positionc - position;
        distf = position - positionf;
        distcc = prod(distc,1);
        distcf = distc(2,:).*distf(1,:);
        distfc = distf(2,:).*distc(1,:);
        distff = prod(distf,1);
        indcc = sub2ind(ssize,positionc(2,:),positionc(1,:));
        indcf = sub2ind(ssize,positionc(2,:),positionf(1,:));
        indfc = sub2ind(ssize,positionf(2,:),positionc(1,:));
        indff = sub2ind(ssize,positionf(2,:),positionf(1,:));
        score = img(indcc).*distcc + img(indcf).*distcf + img(indfc).*distfc + img(indff).*distff;
    end
    score = 1 - (sum(score) / length(points));
    final_score = (((1-frac)*score) + frac*acc);
    if nargout > 1
      parts = [score acc];
    end

end

%***********************************************************
function [position] = get_spline_positions(spline,points,imgsize)
  position = fnval(spline,points);

  indx = (position(1,:) < 1) | (position(1,:) > imgsize(2));
  indy = (position(2,:) < 1) | (position(2,:) > imgsize(1));

  indxy = indx | indy;

  position(:,indxy) = NaN;
end

%***********************************************************
function [surfx,surfy] = get_spline_avg(knots,shift,mymovie,indx)
  polar_knots = carth2elliptic(knots, mymovie.centers(:,indx),mymovie.axes_length(:,indx),mymovie.orientations(:,indx));
  shift = shift ./ sqrt((mymovie.axes_length(1,indx) *cos(polar_knots(:,1))).^2 + (mymovie.axes_length(2,indx) *sin(polar_knots(:,2))).^2);
  inner_pos = elliptic2carth([polar_knots(:,1) polar_knots(:,2)-shift], mymovie.centers(:,indx),mymovie.axes_length(:,indx),mymovie.orientations(:,indx));

  %surfx = [inner_pos(:,1) knots(:,1)];
  %surfy = [inner_pos(:,2) knots(:,2)];
  surfx = inner_pos(:,1);
  surfy = inner_pos(:,2);
  
  return;
end

%***********************************************************
function [values, sampling_pos] = get_spline_profile(img,knots,filter,pts,mymovie,indx)

  values = zeros(length(pts),1);
  shift = floor(min(size(filter))/2);
  convimg = imfilter(img,filter,'symmetric');
  %center = mean(knots, 2);
  polar_knots = carth2elliptic(knots, mymovie.centers(:,indx),mymovie.axes_length(:,indx),mymovie.orientations(:,indx));
  polar_spline = polar_splines(polar_knots);

  shift = shift ./ sqrt((mymovie.axes_length(1,indx) *cos(pts)).^2 + (mymovie.axes_length(2,indx) *sin(pts)).^2);

  %title('profile')
  %t = [0:2*pi/npts:2*pi];
  %t = t(1:end-1);
  r = fnval(polar_spline, pts);
  sampling_pos = elliptic2carth([pts r-shift], mymovie.centers(:,indx),mymovie.axes_length(:,indx),mymovie.orientations(:,indx));
  sampling_pos = round(sampling_pos);
  values = convimg(sub2ind(size(img),sampling_pos(:,2),sampling_pos(:,1)));

  %figure;
  %imshow(img);
  %hold on
  %plot(sampling_pos(:,1),sampling_pos(:,2),'-ob');
  %plot(knots(:,1),knots(:,2),'-or');
  %hold off

  %if(nargout>1)
  %  sampling_pos = round(original);
  %  plot_values = img(sub2ind(size(img),sampling_pos(2,:),sampling_pos(1,:)));
  %end

  return

  %nlayers = size(img,3);
  %values = zeros(length(points),nlayers);
  %for layer=1:nlayers
    xvals = zeros(1,length(points));
    tmpvals = zeros(1,length(points));
    %for i=1:length(points)
    %  position = fnval(spline,points(i));
    %  xvals(end+1) = position(1);
    %  position = round(position);
    %  tmpvals(i) = img(position(2),position(1),layer);
    %end

    position = get_spline_positions(spline,points,size(img));
    index = ~(isnan(position(1,:)));
    position = position(:,index);
    points = points(index);

    %position = fnval(spline,points);
    position = round(position);
    xmin = find(position(1,:) == min(position(1,:)));
    position = [position(:,xmin:end) position(:,1:xmin-1)];
    tmpvals = NaN(1,length(index));
    tmpvals(index) = img(sub2ind(size(img),position(2,:),position(1,:),layer*ones(1,length(points))));
    values(:,layer) = tmpvals;
  %end

  return;
end

%***********************************************************
function [mymovie] = rescale_colors(mymovie)
  
  hon = findobj('Enable','on');
  set(hon,'Enable','off');
  %hwait = waitbar(0,'Rescaling images...','Name','CellCoord Info');
  hwait = findobj('Tag','wait');
  waitbar(0, hwait, 'Rescaling Images...');
  set(hwait,'Visible','on');
  [tmpsize, nframes] = size_data(mymovie.data);
  for i=1:nframes

    %img = double(mymovie(i).data);
    img = load_data(mymovie.data,i);

    pixels = reshape(img);
    vals = hist(pixels, 100);
    figure
    hist(pixels, 100);

    min_val = min(min(min(img)));
    max_val = max(max(max(img)));

    new_img = (img - min_val) ./ (max_val - min_val);

    %mymovie(i).data = new_img;
    mymovie.data = modify_data(mymovie.data,new_img,i);
    waitbar(i/nframes,hwait);
  end
  %close(hwait);
  set(hwait,'Visible','off');
  set(hon,'Enable','on');

  clear img min_val max_val new_img;

  return;
end

%***********************************************************
function [stop] = spline_evol(algo_options, optim_values, flag)

  %fitmov = algo_options.FunctionParameters.fitmov;
  record = false;

  vals = optim_values.x;
  neg_indx = find(vals(1,:)<0);
  vals = abs(vals);
  x = [vals(1,:) vals(1,1)];
  y = [vals(2,:) vals(2,1)];

  [ t,ier] = arcl2d(x,y);
  spline = mycsape(t,[x;y],'periodic');
  [points,distance] = fnplt(spline);
  hwindows = findobj('Tag','evol');
  hax = findobj(hwindows,'Tag','axes');
  hsplines = findobj(hwindows,'LineStyle','-','Marker','none');
  hcp = findobj(hwindows,'LineStyle','none','Marker','o');
  htext = findobj(hwindows,'Tag','text');
  hcurrent = findobj(hwindows,'LineStyle','none','Marker','x');
  if(record)
    hfile = findobj(hwindows,'Tag','filename');
    filename = get(hfile,'String');
  end
  set(hsplines,'XData',points(1,:), 'YData',points(2,:));

  if(strcmpi(flag,'init'))
    line('XData',x(1:end-1),'YData',y(1:end-1),'LineStyle', 'none', 'Marker', 'o', 'Parent',hax, 'Color', 'r');
    set(hwindows,'Visible','on');
    line('XData',x(1),'YData',y(1),'LineStyle', 'none', 'Marker', 'none', 'MarkerSize', 12, 'Parent',hax, 'Color', 'g');
    text(0,0,' ','BackgroundColor',[1 1 1],'Tag','text');
  else
    set(hcp,'XData',x(1:end-1),'YData',y(1:end-1));
    set(htext,'string',['score: ' num2str(optim_values.fval) ', temp: ' num2str(optim_values.temperature(neg_indx))])
    if(~isempty(neg_indx))
      if(isempty(hcurrent))
        hcurrent = findobj(hwindows,'LineStyle','none','Marker','none');
        set(hcurrent,'XData',x(neg_indx),'YData',y(neg_indx),'Marker', 'x');
      else
        set(hcurrent,'XData',x(neg_indx),'YData',y(neg_indx));
      end
    end
  end

  drawnow;
  refresh(hwindows);

  if(record)
    frame = getframe(hax);
    filename = store_data(frame.cdata,filename);
    set(hfile,'String',filename);
  end
  %size(frame.cdata)
  %frame.colormap
  %addframe(fitmov,frame);
  %[img,map] = frame2img(frame);
  %imwrite(img,map,filename,'tiff','WriteMode','append');

  if(strcmpi(flag,'done'))
    set(hwindows,'Visible','off');
    delete(hcp);
    delete(hcurrent);
    delete(htext);
  end
  stop = false;
end

function [color] = find_color(mymovie,img_path)

  f = fields(mymovie);
  for i=1:length(f)
    if(strcmp(mymovie.(f{i}),img_path))
      switch f{i}
        case 'dic'
          color = 'g';
        case 'c100'
          color = 'r';
        case 'c010'
          color = 'b';
        case 'c001'
          color = 'y';
      end
    end
  end
end

%***********************************************************
function [mask] = create_mask(spline,img)

  hon = findobj('Enable','on');
  set(hon,'Enable','off');
  %hwait = waitbar(0,'Creating required mask...','Name','CellCoord Info');
  hwait = findobj('Tag','wait');
  waitbar(0, hwait, 'Creating required mask...');
  set(hwait,'Visible','on');
  [points,dist] = fnplt(spline);

  mask = roipoly(img,points(1,:),points(2,:));
  seround = strel('diamond', 3);
  mask = imerode(mask, seround);

  %close(hwait);
  set(hwait,'Visible','off');
  set(hon,'Enable','on');

  return;
end

function [points,t] = warp_spline(splines, ref_spline, spline_indx, frame_indx, mymovie)

  egg = [splines(1, frame_indx).x' splines(1, frame_indx).y'];
  ref = [ref_spline.x' ref_spline.y'];
  spline = [splines(spline_indx, frame_indx).x' splines(spline_indx, frame_indx).y'];

  %center = mean(egg,2);
  %ref_center = mean(ref,2);
  %figure
  %hold on
  %plot(egg(1,:), egg(2,:),'r');
  %plot(ref(1,:)*10, ref(2,:)*10,'b');
  %plot(spline(1,:), spline(2,:),'g');

  %center = mean(egg,2);
  %ref_center = mean(ref,2);

  [polar_egg] = carth2elliptic(egg, mymovie.centers(:,frame_indx),mymovie.axes_length(:,frame_indx),mymovie.orientations(:,frame_indx));
  [polar_spline] = carth2elliptic(spline, mymovie.centers(:,frame_indx),mymovie.axes_length(:,frame_indx),mymovie.orientations(:,frame_indx));
  [polar_ref] = carth2elliptic(ref, [0; 0],[25; 15],0);

  %figure;
  %plot(polar_egg(:,1),polar_egg(:,2),'-or')
  %hold on;
  %plot(polar_spline(:,1),polar_spline(:,2),'-og')
  %plot(polar_ref(:,1),polar_ref(:,2),'-ob')

  egg = polar_splines(polar_egg');
  %title('egg')
  spline = polar_splines(polar_spline');
  %title('cortex')
  ref = polar_splines(polar_ref');
  %title('ref')

  %pts = polar_spline
  [pts] = fnplt(spline);
  %pts(1,:) = pts(1,:) + (pts(1,:)<0) *2*pi;
  %indx = find(pts(1,:)==min(pts(1,:)));
  %pts = [pts(:,indx:end) pts(:,1:indx-1)];
  %indx = find(pts(1,:)>2*pi,1);
  %if(~isempty(indx))
  %  pts = pts(:,1:indx-1);
  %end
  t = pts(1,:);
  r = pts(2,:);

  valid = (t<=2*pi & t>=0);
  t = t(valid);
  r = r(valid);

  egg_r = fnval(egg, t);
  ref_r = fnval(ref, t);

  warp = r .* ref_r ./ egg_r;
  %plot(t,warp,'-ok')

  if(isfield(mymovie,'isinverted') && mymovie.isinverted)
    points = elliptic2carth([t' warp'], [0;0], [25;15],pi);
  else
    points = elliptic2carth([t' warp'], [0;0], [25;15],0);
  end
  %plot(points(1,:)*10, points(2,:)*10,'k');
  %points = polar2carth([t; ref_r], ref_center);
  %plot(points(1,:)*10, points(2,:)*10,'c');
  %points = polar2carth([t; egg_r], center);
  %plot(points(1,:), points(2,:),'m');
  %if(nargout>1)
  %  original = elliptic2carth([t' r'],  mymovie.centers(:,frame_indx),mymovie.axes_length(:,frame_indx),mymovie.orientations(:,frame_indx));
  %end
  %plot(points(1,:), points(2,:),'y');
end

function [carth_pts] = polar2carth(pts, center)

  [s1 s2] = size(pts);

  if(s1==1)
    pts = pts';
    [s1 s2] = size(pts);
  end
 
  O = pts(1,:);
  r = pts(2,:);

  carth_pts = zeros([s1 s2]);
  carth_pts(1,:,:) = r.*cos(O);
  carth_pts(2,:,:) = r.*sin(O);

  carth_pts = carth_pts + repmat(center, [1 s2]);

  return;
end

function [polar_pts, center] = carth2polar(pts, center)

  [s1 s2] = size(pts);
  if(s2<100)
    %spline = carth_spline(pts(:,[1:end 1]));
    [t,ierr] = arcl2d(pts(1,[1:end 1]),pts(2,[1:end 1]));
    spline = mycsape(t, pts(:,[1:end 1]), 'periodic');

    [new_pts, t] = fnplt(spline);
    pts = fnval(spline, 0:t(end)/100:t(end));
    pts = pts(:,1:end-1);
    [s1 s2] = size(pts);

    center = mean(pts, 2);
  end

  pts = pts - repmat(center, [1 s2]);

  x = squeeze(pts(1,:));
  y = squeeze(pts(2,:));

  polar_pts = zeros([s1 s2]);
  polar_pts(1,:) = atan2(y, x);
  polar_pts(2,:) = sqrt((x.^2)+(y.^2));

  polar_pts(1,:) = polar_pts(1,:) + (2*pi* (polar_pts(1,:)<0));
  indx = find(polar_pts(1,:) == min(polar_pts(1,:)));

  if(indx>1)
    polar_pts = cat(2, polar_pts(:,indx:end), polar_pts(:,1:indx-1));
  end
  if(abs(polar_pts(1,1)-polar_pts(1,end))<pi)
    polar_pts = polar_pts(:,[1 end:-1:2]);
  end

  polar_pts = polar_pts(:,[1:end 1]);
  polar_pts(1,end) = polar_pts(1, end) + 2*pi;

  return;
end

function [splines] = polar_splines(pts)

  if(size(pts,2)==2)
    pts = pts';
  end

  pos = pts(1,:);
  pts = pts(2,:);

  if(size(pos,1)==1)
    pos = pos';
  end
  if(size(pts,1)==1)
    pts = pts';
  end

  twopi = 2*pi;

    if(size(pos,2)==1)
      frame_pos = pos;
    else
      frame_pos = pos(:,i);
    end

    final_pos = frame_pos(1);
    %if(final_pos<0)
      final_pos = final_pos + twopi;
    %else
    %  final_pos = final_pos - twopi;
    %end

    splines = mycsape([frame_pos; final_pos], pts([1:end 1]), 'periodic'); 
    %splines = mycsape([frame_pos(end) - twopi;frame_pos; frame_pos(1) + twopi], pts([end 1:end 1])); 

  return;
end

%***********************************************************
function [x,y,t] = optimize_spline(x, y, img, frac, thresh)

  %x = spline.x;
  %y = spline.y;

  redundant = 1;
  while (redundant>0)
     [t,redundant] = arcl2d(x,y);
     if  redundant > 0
        delpt(redundant); 
     end
  end

  position = [x;y];

  orig_score = get_spline_score([x(1:end-1);y(1:end-1)],img,frac);
  %numpts = length(x);
  numpts = 200;

  new_score = orig_score;
  while (abs(new_score-orig_score)<thresh) && (numpts > 11)
    new_numpts = round(numpts*0.9)-1;

    x = position(1,:);
    y = position(2,:);

    [t,ier] = arcl2d(x,y);
    spline = csape(t,position,'periodic');
    [points,distance] = fnplt(spline);
    spline.dist=distance;

    points = [0:spline.dist(end)/new_numpts:spline.dist(end)];
    
    position = get_spline_positions(spline,points(1:end-1),size(img));
    index = ~(isnan(position(1,:)));
    position = position(:,index);

    xmin = find(position(1,:) == min(position(1,:)));
    position = [position(:,xmin:end) position(:,1:xmin-1)];
  
    new_score = get_spline_score(position(:,1:end),img,frac);

    numpts = new_numpts;
    position = [position position(:,1)];
  end
  return;

    %position = get_spline_positions(spline,points,size(img));
    %index = ~(isnan(position(1,:)));
    %position = position(:,index);
    %points = points(index);

    %position = fnval(spline,points);
    %position = round(position);
    %xmin = find(position(1,:) == min(position(1,:)));
    %position = [position(:,xmin:end) position(:,1:xmin-1)];
    %tmpvals = NaN(1,length(index));
    %tmpvals(index) = img(sub2ind(size(img),position(2,:),position(1,:),layer*ones(1,length(points))));
    %values(:,layer) = tmpvals;

          %[ttmp,ier] = arcl2d([bestx bestx(1)],[besty besty(1)]);
          %spline = csape(ttmp,[bestx bestx(1);besty besty(1)],'periodic');
          %[points,distance] = fnplt(spline);
          %tmpscore = get_spline_score(img,spline,[0:distance(end)/1000:distance(end)],'direct','subpixsum');

  dist = diff(t);
  while (any(dist < min_dist))
    [ordered_dist,index] = sort(dist);
    delpt(index(1));
    [t,redundant] = arcl2d(x,y);
    dist = diff(t);
  end

  return;

function delpt(i)

    j = i+1:length(x);               % indices of control points
    x(j-1) = x(j); 
    y(j-1) = y(j);
    x(end) = [];
    y(end) = [];

    if i == 1
       x(end) = x(1);
       y(end) = y(1);
    end

    return;
  end
end
