function mymovie = shapes2movie(mymovie, shape)

  isspline = isfield(shape,'form');
  [imgsize, nframes] = size_data(mymovie.c111);

  for i=1:nframes-1
    if (i<size(shape,2))
      if (~isspline)
        mymovie.egg(i).path = shape(1,i).path;
        mymovie.cortex(i).path = shape(2,i).path;
        
        mymovie.egg(i).spline = create_spline(mymovie.egg(i).path);
        mymovie.cortex(i).spline = create_spline(mymovie.cortex(i).path);
      else
        mymovie.egg(i).spline = shape(1,i);
        mymovie.cortex(i).spline = shape(2,i);
      end
        
      [center, axes_length, orientation] = fit_ellipse(mymovie.egg(i).spline);
      %mymovie.egg(i).spline = realign_spline(mymovie.egg(i).spline, center, orientation);
      %mymovie.cortex(i).spline = realign_spline(mymovie.cortex(i).spline, center, orientation);

      mymovie.centers(:,i) = center';
      mymovie.axes_length(:,i) = axes_length';
      mymovie.orientations(i) = orientation;
    else
      mymovie.egg(i).path = [];
      mymovie.cortex(i).path = [];
    end
  end

