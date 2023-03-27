function centrosome_movie(mymovie, fname)

  [imgsize, nframes] = size_data(mymovie.data);
  himg = [];
  hcentr1 = [];
  hcentr2 = [];

  for i=1:nframes
    nimg = i;

    img = imnorm(double(load_data(mymovie.data, nimg)));
    img = realign(img, [480 640], mymovie.data.centers(:,nimg), mymovie.data.orientations(1,nimg));
    centr1 = realign(mymovie.data.centrosomes(nimg).carth(1,:),[480 640],mymovie.data.centers(:,nimg), mymovie.data.orientations(1,nimg)).';
    centr2 = realign(mymovie.data.centrosomes(nimg).carth(2,:),[480 640],mymovie.data.centers(:,nimg), mymovie.data.orientations(1,nimg)).';

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
