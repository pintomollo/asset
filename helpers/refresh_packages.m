function refresh_packages()

  pkg unload image;
  pkg unload statistics;
  pkg unload control;
  pkg unload matgeom;
  pkg unload struct;
  pkg unload optim;
  pkg unload signal;
  pkg unload tisean;
  pkg unload io;

  pkg load io;
  pkg load statistics;
  pkg load image;

  return;
end
