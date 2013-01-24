function [signals, full_pos, shuffle] = combine_domains(mymovies, sync_type, sync_param, min_counts)

  if (nargin == 1)
    sync_type = 'lsr';
    sync_param = [];
    min_counts = 0.2;
  elseif (nargin == 2 & isnumeric(sync_type))
    min_counts = sync_type;
    sync_type = 'lsr';
    sync_param = [];
  end

  if (ischar(mymovies))
    if (~isempty(findstr(mymovies, '*')))
      tmp = dir(mymovies); 
      mymovies = cell(1, length(tmp));

      for i=1:length(tmp)
        mymovies{i} = tmp(i).name;
      end
    else
      mymovies = {mymovies}; 
    end
  elseif (~iscell(mymovies))
    tmp = mymovies;
    mymovies = {};
    mymovies{1} = tmp;
  end

  nmovies = length(mymovies);

  signals = [];
  align_time = 0;

  shuffle = randperm(nmovies);
  for i=shuffle
    load(mymovies{i});
    
    [fraction, max_width, cell_width, domain, pos] = domain_expansion(mymovie, opts);
    time = get_manual_timing(mymovie, opts);
    domain = domain(1:time(end), :).';

    size_diff = size(signals, 1) - size(domain, 1);
    if (size_diff > 0)
      domain = padarray(domain, [size_diff/2 0], NaN);
    end
    if (~isempty(signals) & size_diff < 0)
      signals = padarray(signals, [-size_diff/2 0 0], NaN);
    end

    switch (sync_type)
      case 'fraction'
        align_time = find(fraction >  sync_param, 1, 'first');
        if (isempty(align_time))
          align_time = time(end);
        end
      case 'lsr'
        domain = domain / mymean(domain(:));
        [signals, new_time] = find_min_residue(signals, align_time, domain, time(1), 0.75, 20);
        if (new_time == 0)
          keyboard
        end
        align_time = new_time;
    end
  end

  keyboard

  [h,w,n] = size(signals);
  full_pos = [1:w] - (((w-1)/2)+1);
  full_pos = full_pos * (median(diff(pos)));

  counts = sum(~isnan(signals), 3);
  if (min_counts < 1)
    min_counts = round(min_counts*n);
  end
  rows = any(counts >= min_counts, 2);
  cols = any(counts >= min_counts, 1);

  signals = signals(rows, cols, :);
  full_pos = full_pos(cols);

  return;
end
