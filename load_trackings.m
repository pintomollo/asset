function [trackings, max_frames] = load_trackings(trackings, opts)

  max_frames = 0;
  if (isfield(trackings, 'child'))
    for i=1:length(trackings.child)
      [trackings.child(i), nframes] = load_trackings(trackings.child(i), opts);
      if (nframes > max_frames)
        max_frames = nframes;
      end
    end
  end

  if (isfield(trackings, 'expr'))
    trackings.expr = relativepath(trackings.expr);
  end

  for i=1:length(trackings.files)
    if (isempty(trackings.files(i).shapes) || opts.recompute)
      trackings.files(i).fname = relativepath(trackings.files(i).fname);

      [trackings.files(i).shapes, trackings.files(i).groups] = load_shapes(trackings.files(i).fname);
      trackings.files(i).splines = shapes2splines(trackings.files(i).shapes);
    end
    if (opts.recompute)

      if (opts.follow_periphery)
        trackings.files(i).splines = spline_periphery(trackings.files(i).splines);
      end

      nframes = size(trackings.files(i).shapes,2);
      if (nframes > max_frames)
        max_frames = nframes;
      end
    end
  end

  return;
end
