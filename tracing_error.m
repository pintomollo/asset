function mymovie = tracing_error(trackings, mymovie, field, opts)

  mymovie.(field).splines = movie2spline(mymovie, field, opts);

  if (opts.follow_periphery)
    mymovie.(field).splines = spline_periphery(mymovie.(field), opts);
  end

  mymovie.(field).errors = spline_error(trackings.(field), mymovie.(field), opts.nbins, opts);
  mymovie.(field).update = false(size(mymovie.(field).update));

  return;
end
