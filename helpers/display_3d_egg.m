function display_3d_egg(mymovie, opts, save_frames)

  if (nargin == 1)
    opts = get_struct('ASSET');
    save_frames = false;
  elseif (nargin == 2)
    if (islogical(opts))
      save_frames = opts;
      opts = get_struct('ASSET');
    else
      save_frames = false;
    end
  end

  if (ischar(mymovie))

    fname = mymovie;
    tmp = dir(fname); 

    for i=1:length(tmp)
      name = tmp(i).name
      load(name);

      display_3d_egg(mymovie, opts, save_frames);
    end
    
    return;
  end

  if (save_frames)
    indxs = findstr(mymovie.experiment(1:end-1), '_');
    movie_number = mymovie.experiment(indxs(end)+1:end-1);
  end
  
  [nframes, imgsize] = size_data(mymovie.dic);
  for i=1:nframes
    nimg = i;
    
    z_coef = sqrt(1 - ((mymovie.metadata.z_position(i) - mymovie.metadata.center_3d(3)).^2 / mymovie.metadata.axes_length_3d(3).^2));

    hold off;
    imshow(imnorm(double(load_data(mymovie.dic, nimg))));
    hold on;
    draw_ellipse(mymovie.dic.centers(:, nimg), mymovie.dic.axes_length(:, nimg), mymovie.dic.orientations(1, nimg));
    if (isreal(z_coef))
      draw_ellipse(mymovie.metadata.center_3d(1:2) / opts.pixel_size, mymovie.metadata.axes_length_3d(1:2) * z_coef / opts.pixel_size, mymovie.metadata.orientation_3d, 'r');
    end

    if (save_frames)
      print('-dpng', ['./PNG/Z-size-' movie_number  '-' num2str(i) '.png']);
    else
      pause;
    end
  end


  return;
end
