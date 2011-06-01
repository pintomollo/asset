function mymovie = identify_channels(channels, opts)
% IDENTIFY_CHANNELS identifies the type of data contained in each channel of the movie,
% including the type of filtering required and prompts the user for verification. In
% addition, this function incorporates the different channels into the "mymovie"
% structure (see get_struct.m).
%
%   [MYMOVIE] = IDENTIFY_CHANNEL(CHANNELS, OPTS) converts the channels, initially an
%   array of "channel" structures (see get_struct.m), into the mymovie structure, while
%   inferring the type of data in each of them.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 19.05.2011

  % Get the number of channels and loop over each of them to make an initial guess
  nchannels = length(channels);
  for i=1:nchannels
    
    % Basic idea is based on my data: mCherry-PH & GFP-anything_else.
    % Consequently the mCherry channel is supposedly the cortex. 
    if (strfind(channels(i).file, 'mCherry'))
      channels(i).color = [1 0 0];
      channels(i).type = 'cortex';

    % And the GFP the data
    elseif (strfind(channels(i).file, 'FITC') | strfind(channels(i).file, 'GFP'))
      channels(i).color = [0 1 0];
      channels(i).type = 'data';
    end
  end

  % Get THE structure
  mymovie = get_struct('mymovie', 1);

  % In case the verbosity is null, do not prompt for confirmation, useful for batch processing
  if (opts.verbosity > 1)
    channels = input_channels(channels);
  end

  % Copy each channel into mymovie
  for i=1:nchannels

    % Get the type of data as this defines its position in mymovie
    channel_name = channels(i).type;

    % Copy the data
    if (~isfield(mymovie, channel_name) | isempty(mymovie.(channel_name)))
      mymovie.(channel_name) = channels(i);
    else
      mymovie.(channel_name)(end+1) = channels(i);
    end
  end

  return;
end
