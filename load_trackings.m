function [trackings, max_frames] = load_trackings(trackings, opts)


%%%%%%%%%%%%%% DO NOT USE THE MEAN FIELD, USE THE eggshell AND cortex RATHER

  max_frames = 0;
  if (~isfield(trackings, 'child'))
    fields = fieldnames(trackings);

    for i=1:length(fields)
      field = fields{i};
      if (isstruct(trackings.(field)))
        trackings.(field) = load_trackings(trackings.(field), opts, hwait);
      end
    end

    return;
  end

  if (isfield(trackings.child, 'shapes'))
    for i=1:length(trackings.child)
      if (isempty(trackings.child(i).shapes) || opts.recompute)
        trackings.child(i).fname = relativepath(trackings.child(i).fname);

        [trackings.child(i).shapes, trackings.child(i).groups] = load_shapes(trackings.child(i).fname);
      end
      if (opts.recompute)

        if (opts.follow_periphery)
          trackings.child(i).shapes = path_periphery(trackings.child(i).shapes);
        end

        nframes = size(trackings.child(i).shapes,2);
        if (nframes > max_frames)
          max_frames = nframes;
        end
      end
    end
  else
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

  return;
end
