function [mymovie, opts] = filter_channels(mymovie, opts)

  all_channels = {};

  nfilters = length(opts.filters);
  filters = cell(nfilters, 1);
  channels = cell(nfilters, 1);
  params = cell(nfilters, 1);
  boundaries = cell(nfilters, 1);

  boundaries(:) = {'symmetric'};
  new_movie = mymovie;

  for i=1:nfilters
    filter = opts.filters(i);

    if (filter.applied | isempty(filter.filter))
      continue;
    else
      opts.filters(i).applied = true;
    end

    if (~iscell(filter.params))
      filter.params = {filter.params};
    end
    if (~isempty(filter.params) & ischar(filter.params{end}))
      boundaries{i} = filter.params{end};
      filter.params = filter.params(1:end-1);
    end
    if (ischar(filter.filter))
      new_movie.experiment = [new_movie.experiment filter.filter '_'];

      if (exist(filter.filter, 'file') ~= 0)
        func = str2func(filter.filter);

        filters{i} = func;
        params{i} = filter.params;
        if (strncmp(filter.filter, 'binning', 7))
          channels{i} = 'all';
          opts.binning = filter.params{1};
        else
          channels{i} = filter.channel;
        end
      else
        try
          h = fspecial(filter.filter, filter.params{:});
        catch ME
          warning(['Filter ' filter.filter ' is unknown, ignoring it']);

          continue;
        end
        filters{i} = h;
        channels{i} = filter.channel;
      end
    elseif (isa(filter.filter, 'function_handle'))
      filters{i} = filter.filter;
      params{i} = filter.params;
      channels{i} = filter.channel;

      new_movie.experiment = [new_movie.experiment func2str(filter.filter) '_'];
    elseif (isnumeric(filter.filter))
      filters{i} = filter.filter;
      channels{i} = filter.channel;

      new_movie.experiment = [new_movie.experiment 'custom_filter_'];
    else
      warning(['Do not know what to with a filter of type ' class(filter.filter)]);
    end

    if (strncmp(channels{i}, 'all', 3))
      if (isempty(all_channels))
        all_channels = get_channels(mymovie);
      end

      for j=1:length(all_channels)
        if (j==1)
          indx = i;
        else
          indx = length(filters) + 1;
        end

        channels(indx) = all_channels(j);
        filters(indx) = filters(i);
        params(indx) = params(i);
        boundaries(indx) = boundaries(i);
      end
    end
  end

  valids = (~cellfun('isempty',filters));
  filters = filters(valids);
  channels = channels(valids);
  params = params(valids);
  boundaries = boundaries(valids);
  nfilters = length(filters);

  if (nfilters > 0)
    nframes = size_data(mymovie.(channels{1}));

    for i = 1:nframes
      prev_channel = 'none';
      for f = 1:nfilters
        if (~strcmp(prev_channel, channels{f}))
          if (f > 1)
            save_data(new_movie.(prev_channel), img);
          end
          img = load_data(mymovie.(channels{f}), i);
          prev_channel = channels{f};

          if (i == 1)
            new_movie.(prev_channel).fname = absolutepath(get_new_name('tmpmat(\d+)\.ome\.tiff?', 'TmpData'));
          end
        end
        if (isnumeric(filters{f}))
          img = imfilter(img, filters{f}, boundaries{f});
        elseif (isa(filters{f}, 'function_handle'))
          img = filters{f}(img, params{f}{:});
        end
      end
      save_data(new_movie.(prev_channel), img);
    end

    mymovie = new_movie;
  end

  return;
end
