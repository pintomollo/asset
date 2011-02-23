function centrosome_movie(mymovie, fname)

  [imgsize, nframes] = size_data(mymovie.dic);
  himg = [];
  hcentr1 = [];
  hcentr2 = [];

  for i=1:nframes
    nimg = i;

    img = load_data(mymovie.data, i);
    img = realign(img, [480 640], mymovie.dic.centers(:,i), mymovie.dic.orientations(1,i));
    centr1 = realign(mymovie.data.centrosomes(1).position([2:-1:1],i),[480 640],mymovie.dic.centers(:,nimg), mymovie.dic.orientations(1,nimg)).';
    centr2 = realign(mymovie.data.centrosomes(2).position([2:-1:1],i),[480 640],mymovie.dic.centers(:,nimg), mymovie.dic.orientations(1,nimg)).';

    if (isempty(himg))
      himg = imshow(img);
      hold on
      hcentr1 = myplot(centr1, 'LineStyle','none','Marker','o','MarkerFaceColor',[1 0 0]);
      hcentr2 = myplot(centr2, 'LineStyle','none','Marker','o','MarkerFaceColor',[0 0 1]);
    else
      set(himg, 'CData', img);
      set(hcentr1, 'XData', centr1(1), 'YData', centr1(2));
      set(hcentr2, 'XData', centr2(1), 'YData', centr2(2));
    end

    drawnow;
    movie(i) = getframe();
  end

  movie2avi(movie,fname,'FPS',6);

  return;
end
