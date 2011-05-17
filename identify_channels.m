function mymovie = identify_channels(channels, opts)

  ncolors=3;
  nchannels = length(channels);
  detrend = zeros(nchannels,1);

  for i=1:nchannels
    if (strfind(channels(i).file, 'mCherry'))
      channels(i).color = [1 0 0];
      channels(i).type = 'cortex';
    elseif (strfind(channels(i).file, 'FITC'))
      channels(i).color = [0 1 0];
      channels(i).type = 'data';
    end
  end

  mymovie = get_struct('mymovie', 1);
  if (opts.verbosity > 1)
    channels = input_channels(channels);
  end

  for i=1:nchannels
    channel_name = channels(i).type;

    if (~isfield(mymovie, channel_name) | isempty(mymovie.(channel_name)))
      mymovie.(channel_name) = channels(i);
    else
      mymovie.(channel_name)(end+1) = channels(i);
    end
  end
end
