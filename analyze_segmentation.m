function [mymovie, trackings] = analyze_segmentation(mymovie, trackings, opts)

  if (strncmp(opts.segmentation_type, 'all', 3))
    types = {'dic', 'markers'};
  else
    types = {opts.segmentation_type};
  end

  for t=1:length(types)
 
    if (isempty(trackings.(types{t}).average) || opts.recompute)
      trackings.(types{t}) = mean_trackings(trackings.(types{t}), opts);
    end

    if (isempty(trackings.(types{t}).errors) || opts.recompute)
      trackings.(types{t}) = tracking_error(trackings.(types{t}), opts);
    end

    mymovie = tracing_error(trackings, mymovie, types{t}, opts);

    if (opts.auto_save)
      save(mymovie.experiment, 'mymovie', 'trackings');
    end
  end

  return;
end
