function display_quantification(mymovie, fname)

  if (nargin == 1)
    fname = '';
  end

  nframes = length(mymovie.data.quantification);

  if (size(mymovie.data.quantification(1).cortex, 1) == 6 | size(mymovie.data.quantification(1).cortex, 2) == 6)
    quant = 'front';

    max_val = 0;
    for i=1:nframes
      tmp_max = max(mymovie.data.quantification(i).cortex(:,3));

      if (tmp_max > max_val)
        max_val = tmp_max;
      end
    end
  else
    quant = 'line';
  end
  
  pos = [-16:16];

  for i=1:nframes
    switch quant
      case 'line'

    pts = mymovie.data.cortex(i).carth;
    ell_pts = carth2elliptic(pts, mymovie.data.centers(:, i), mymovie.data.axes_length(:,i), mymovie.data.orientations(1,i), true);

    if (i == 1)
      hfig = figure;

      hax = axes;
    set(hax, ...
    'FontSize',18, ...
    'NextPlot','add', ...
    'Tag', 'axes',...
    'Parent',hfig, ...
    'Xlim', [0 2*pi]);


      signal = line(ell_pts(:,1),mymovie.data.quantification(i).cortex, ...
                'Color', [0 0.5 0], ...
                'HandleVisibility', 'callback', ...
                'EraseMode', 'none', ...
                'HitTest','off', ...
                'LineStyle','-', 'Marker','none', ...
                'LineWidth', 2, ...
                'Parent',hax, ...
                'SelectionHighlight','off', ...
                'Tag','Eggshell');


    else
      set(signal, 'XData', ell_pts(:,1), 'YData', mymovie.data.quantification(i).cortex);
    end

        refresh(hfig);
        drawnow;
        if (~isempty(fname))
          movie(i) = getframe(hfig);
        end

    case 'front'
      params = mymovie.data.quantification(i).cortex;
      tmp_params = params;
      tmp_params(:, [4 end]) = 0;
 
      if (i == 1)
        figure;
      end
      subplot(221)
      imagesc(front_function(tmp_params, pos), [0 max_val]);

      tmp_params = params;
      tmp_params(:, [3 4]) = 0;

      subplot(222)
      imagesc(front_function(tmp_params, pos), [0 max_val]);

      subplot(223)
      imagesc(front_function(params, pos), [0 max_val]);
      title(num2str(i));

      pause
    end
  end

  if (~isempty(fname))
    movie2avi(movie,fname,'FPS',15);
  end

  return;
end
