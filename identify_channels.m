function mymovie = identify_channels(channels)

  ncolors=3;
  nchannels = length(channels);
  detrend = zeros(nchannels,1);

  mymovie = get_struct('mymovie', 1);
  channels = input_channels(channels);

  for i=1:nchannels
    channel_name = channels(i).type;

    if (~isfield(mymovie, channel_name) | isempty(mymovie.(channel_name)))
      mymovie.(channel_name) = channels(i);
    else
      mymovie.(channel_name)(end+1) = channels(i);
    end
  end
end
