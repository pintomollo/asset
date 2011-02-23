function display_movie(mymovie, opts, fname)

  if (nargin < 3)
    fname = '';
  end

  crop = false;
  type = opts.segmentation_type;
  cyto = mymovie.(type).cytokinesis;
  if (opts.crop_export)
    crop = max(mymovie.(type).axes_length, [], 2) * opts.crop_size;
    crop = round(crop([2 1],1).');
  end

  fields = {'dic', 'eggshell', 'cortex'; ...
            'img', 'carth',    'carth'};
  switch type
    case 'markers'
      if (isfield(mymovie, 'cortex') & ~isempty(mymovie.cortex))
        fields{1,1} = 'cortex';
      else
        fields{1,1} = 'dic';
      end
      fields = [fields {'ruffles' 'ruffles'; 'carth' 'paths'}];
    case 'data'
      fields{1,1} = 'data';
      fields = [fields {'centrosomes'; 'carth'}];
  end

  [imgsize, nframes] = size_data(mymovie.(fields{1,1}));
  if (all(crop))
    imgsize = crop;
  end
  handles = struct('img', -1, ...
                   'eggshell', -1, ...
                   'cortex', -1, ...
                   'ruffles', -1, ...
                   'centrosomes', [-1 -1]);

  mygray = [0:255]' / 255;
  mygray = [mygray mygray mygray];

  hfig = figure('Colormap',mygray);
  haxes = axes('Parent', hfig, 'DataAspectRatioMode', 'manual', 'DataAspectRatio', [1 1 1]);
  fps = 6;
  t = timer('Period', 1/fps, 'TimerFcn', @(varargin){});
  start(t);

  for nimg = 1:nframes
    values = extract_fields(mymovie, nimg, type, fields, crop);
    curr_time = (nimg - cyto) * 10;
    curr_time = [floor(curr_time/60) abs(rem(curr_time, 60))];
    if (curr_time(2) == 0)
      if (curr_time(1) < 0)
        curr_time(1) = curr_time(1) - 1;
      end
      curr_time = [num2str(curr_time(1)) ':00'];
    else
      curr_time = [num2str(curr_time(1)) ':' num2str(curr_time(2),2)];
    end

    if (nimg == 1)
      handles.img = image(values{1},'Parent',haxes,'CDataMapping','scaled');
      set(haxes,'Visible','off');
      handles.eggshell = line(values{2}(:,1),values{2}(:,2),'Color',[0 1 0],'Parent',haxes);
      handles.cortex = line(values{3}(:,1),values{3}(:,2),'Color',[1 0.5 0],'Parent',haxes);
      switch type
        case 'markers'
          handles.ruffles(1) = line(values{4}(:,1),values{4}(:,2),'Marker','d','Color',[1 0 1],'LineStyle','none','Parent',haxes, 'MarkerFaceColor', [1 0 1]);
          handles.ruffles(2) = line(values{5}(:,1),values{5}(:,2),'Parent',haxes, 'Color', [1 0.5 0]);
        case 'data'
          handles.centrosomes(1) = line(values{4}(1,1),values{4}(1,2),'Marker','o','Color',[0 0 1],'Parent',haxes);
          handles.centrosomes(2) = line(values{4}(2,1),values{4}(2,2),'Marker','o','Color',[1 0 0],'Parent',haxes);
      end

      handles.time = text(imgsize(2),imgsize(1),0,curr_time,'BackgroundColor',[1 1 1],'VerticalAlignment','bottom', 'HorizontalAlignment','right','FontSize', 20,'Margin',1, 'Parent', haxes);

      set(haxes, 'DataAspectRatio', [1 1 1]);
    else
      set(handles.img, 'CData', values{1});
      set(handles.eggshell,'XData',values{2}(:,1),'YData',values{2}(:,2));
      set(handles.cortex,'XData',values{3}(:,1),'YData',values{3}(:,2));

      switch type
        case 'markers'
          set(handles.ruffles(1), 'XData', values{4}(:,1), 'YData', values{4}(:,2));
          set(handles.ruffles(2), 'XData', values{5}(:,1), 'YData', values{5}(:,2));
        case 'data'
          set(handles.centrosomes(1), 'XData', values{4}(1,1), 'YData', values{4}(1,2));
          set(handles.centrosomes(2), 'XData', values{4}(2,1), 'YData', values{4}(2,2));
      end

      set(handles.time, 'String', curr_time);
    end

    drawnow;

    if (~isempty(fname))
      movie(nimg) = getframe(haxes);
    else
      wait(t);
      start(t);
    end
  end

  if (~isempty(fname))
    movie2avi(movie,fname,'FPS',fps);
  end

  return;
end
