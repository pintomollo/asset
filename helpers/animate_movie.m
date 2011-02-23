function animate_movie(mymovie, trackings, fname, opts, displ_type)

  if (nargin < 4)
    opts = get_struct('RECOS',1);
  end
  
  type = opts.segmentation_type;

  if (nargin < 5)
    displ_type = type;
  end

  [imgsize, nframes] = size_data(mymovie.(displ_type));
  hline = [];
  hpts = [];
  hcentr = [];

  for i=1:nframes
    img = load_data(mymovie.(displ_type), i);
    egg = mymovie.(type).eggshell(i).carth;
    cortex = mymovie.(type).cortex(i).carth;
    if (~isempty(trackings))
      megg = trackings.(type).mean(1,i);
      mcortex = trackings.(type).mean(2,i);
    end
    if (isfield(mymovie.(type), 'ruffles'))
      pts = [mymovie.(type).ruffles(i).carth; NaN(1,2)];
    else
      pts = [];
    end
    if (isfield(mymovie, 'data') & isfield(mymovie.data, 'centrosomes'))
      centr = [mymovie.data.centrosomes(i).carth; zeros(1,2)];
    else
      centr = [];
    end


    if (opts.crop_export)
      [s,c,o] = deal([388 591], mymovie.(type).centers(:,i), mymovie.(type).orientations(1,i));
      img = realign(img, s, c, o);
      egg = realign(egg, s, c, o);
      cortex = realign(cortex, s, c, o);

      if (~isempty(trackings))
        megg = realign(megg, s, c , o);
        mcortex = realign(mcortex, s, c, o);
      end
      if (~isempty(pts))
        pts = realign(pts, s, c, o);
      end
      if (~isempty(centr))
        centr = realign(centr, s, c, o);
      end
    end
    
    if (i == 1)
      himg = imshow(img);
      hold on;

      if (~isempty(trackings))
        h = myplot(megg, 'b', 'LineWidth', 2);
        hline = [hline h];
        h = myplot(mcortex, 'r', 'LineWidth', 2);
        hline = [hline h];
      end

      h = myplot(egg, 'Color', [0 1 0], 'LineWidth', 2);
      hline = [hline h];
      h = myplot(cortex, 'Color', [1 0.5 0], 'LineWidth', 2);
      hline = [hline h];

      if (~isempty(pts))
        hpts = myplot(pts, 'Color', [1 0 0], 'Marker', 'o', 'LineStyle', 'none', 'LineWidth', 2);
      end
      if (~isempty(centr))
        h = myplot(centr(:,1), 'Color', [1 0 0], 'Marker', 'o', 'LineStyle', 'none', 'LineWidth', 2);
        hcentr = [hcentr h];
        h = myplot(centr(:,2), 'Color', [0 0 1], 'Marker', 'o', 'LineStyle', 'none', 'LineWidth',2);
        hcentr = [hcentr h];
      end
    else
      set(himg, 'CData', img);

      indx = 1;
      if (~isempty(trackings))
        h = myplot(megg, hline(indx:indx+1));
        indx = indx + length(h);
        h = myplot(mcortex, hline(indx:indx+1));
        indx = indx + length(h);
      end
      h = myplot(egg, hline(indx:indx+1));
      indx = indx + length(h);
      h = myplot(cortex, hline(indx:end));

      if (~isempty(pts))
        h = myplot(pts, hpts);
      end

      if (~isempty(centr))
        h = myplot(centr(1,:), hcentr(1));
        h = myplot(centr(2,:), hcentr(2));
      end
    end

    drawnow;
    if (~isempty(fname))
      movie(i) = getframe();
    end
  end

  if (~isempty(fname))
    movie2avi(movie,fname,'FPS',6);
  end

  return;
end
