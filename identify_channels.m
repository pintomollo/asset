function mymovie = identify_channels(channels)

  ncolors=3;
  nchannels = length(channels);
  detrend = zeros(nchannels,1);

  mymovie = get_struct('mymovie', 1);

  bw = true;
  update = false;

  maxframes=Inf;
  for i=1:nchannels
    [imgsize nframes] = size_data(channels(i));
    if (nframes<maxframes)
      if (~isinf(maxframes))
        update = true;
      end
      maxframes=nframes;
    end

    maxchan = 0;
    colors=zeros(1,ncolors);
    for j=1:ncolors
      img = load_data(channels(i), j);
      colors(j) = sum(sum(img>0))/prod(imgsize);
    end
    channels(i).color = double((max(colors)-colors)<0.01);

    if (all(channels(i).color))
      channels(i).type = 'dic';
    elseif (all(channels(i).color == [1 0 0]))
      channels(i).type = 'cortex';
      channels(i).detrend = true;
      bw = false;
    elseif (all(channels(i).color == [0 0 1]))
      channels(i).type = 'eggshell';
      channels(i).detrend = true;
      bw = false;
    else
      bw = false;
    end
  end

  channels = input_channels(channels);
  data_ctr = 0;

  for i=1:nchannels
    img=zeros(imgsize);

    if (strcmp(channels(i).type, 'data'))
      data_ctr = data_ctr + 1;
      indx = data_ctr;
    else
      indx = 1;
    end
    channel_name = channels(i).type;

    if (indx == 1)
      mymovie.(channel_name) = channels(i);
    else
      mymovie.(channel_name)(indx) = channels(i);
    end
    if (~bw | update)

      mymovie.(channel_name)(indx).fname = '';
      for j=1:maxframes
        if (bw)
          mymovie.(channel_name)(indx) = store_data(mymovie.(channel_name)(indx), load_data(channels(i),j));
        else

          if(mod(j,ncolors)~=0)
            img = img + load_data(channels(i), j);
          else
            mymovie.(channel_name)(indx) = store_data(mymovie.(channel_name)(indx), img / ncolors);
            img = zeros(imgsize);
          end
        end
      end

      delete(channels(i).fname);
    end
  end
end
