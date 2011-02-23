function mymovie = perform_application(mymovie, opts)

  if (ischar(opts.application))
    opts.application = {opts.application};
  end

  for i=1:length(opts.application)
    if (~isempty(opts.application{i}))
      switch opts.application{i}
        case {'ruffling', 'ruffles'}
          mymovie = track_ruffles(mymovie, opts);
        case 'centrosomes'
          mymovie = track_centrosomes(mymovie, opts);
        case 'quantification'
          mymovie = quantify_signal(mymovie, opts);
        otherwise
          warning(['Application ' opts.application{i} ' not implemented.']);
      end
    end
  end

  return;
end
