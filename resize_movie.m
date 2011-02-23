function mymovie = resize_movie(old_movie, framerate, imgscale)

  [imgsize, nframe] = size_data(old_movie.c111);
  mymovie = old_movie;

  if (nargin < 3)
    imgscale = 1;
  end

  mymovie.c111 = '';
  mymovie.c100 = '';
  mymovie.c010 = '';
  mymovie.c001 = '';
  if (isfield(mymovie,'polar'))
    mymovie.polar = '';
  end

  for i=1:framerate:nframe
    mymovie.c111 = store_data(imresize(load_data(old_movie.c111,i),imgscale), mymovie.c111);
    mymovie.c100 = store_data(imresize(load_data(old_movie.c100,i),imgscale), mymovie.c100);
    mymovie.c010 = store_data(imresize(load_data(old_movie.c010,i),imgscale), mymovie.c010);
    mymovie.c001 = store_data(imresize(load_data(old_movie.c001,i),imgscale), mymovie.c001);

    if (isfield(mymovie,'polar'))
      mymovie.polar = store_data(imresize(load_data(old_movie.polar,i),imgscale), mymovie.polar);
    end
  end

  return;
end
