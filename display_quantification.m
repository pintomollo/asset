function display_quantification(mymovie, fname)


  if (nargin == 1)
    fname = '';
  end

  nframes = length(mymovie.data.quantification);

  for i=1:nframes
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
    'Ylim', [0 1], ...
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
  end

  if (~isempty(fname))
    movie2avi(movie,fname,'FPS',15);
  end

  return;
end
