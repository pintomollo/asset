function mymovie = tracing_error(trackings, mymovie, field, opts)

  % Do we need to store this ?
  segments = extract_segmentations(mymovie, field, opts);
  update = mymovie.(field).update;

  if (opts.follow_periphery)
    segments = path_periphery(segments, update, opts);
  end

  if (isfield(mymovie.(field), 'errors') & ~isempty(mymovie.(field).errors))
    mymovie.(field).errors = path_error(trackings.(field).average, segments, mymovie.(field).errors, mymovie.(field).update, trackings.(field).reference, opts);
  else
    mymovie.(field).errors = path_error(trackings.(field).average, segments, mymovie.(field).update, trackings.(field).reference, opts);
  end
  mymovie.(field).update = false(size(mymovie.(field).update));

  return;
end
