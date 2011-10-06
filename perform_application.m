function mymovie = perform_application(mymovie, opts)

  if (ischar(opts.application))
    opts.application = {opts.application};
  end

  done_timing = false;

  for i=1:length(opts.application)
    if (~isempty(opts.application{i}))
      switch opts.application{i}
        case {'ruffling', 'ruffles'}
          if (~done_timing)
            mymovie = track_ruffles(mymovie, opts);
          end
        case 'centrosomes'
          mymovie = track_centrosomes(mymovie, opts);
        case 'quantification'
          mymovie = cortical_signal(mymovie, opts);
        case 'timing'
          mymovie = time_cell_cycle(mymovie, opts);
          done_timing = true;
        case 'z-size'
          mymovie = reconstruct_egg(mymovie, opts);
        case 'nuclei'
          if (~done_timing)
            mymovie = detect_dic_nuclei(mymovie, opts);
          end
        otherwise
          warning(['Application ' opts.application{i} ' not implemented.']);
      end
    end
  end

  return;
end
