function [signals, full_pos, shuffle] = combine_domains(mymovies, sync_type, sync_param, min_counts)

  %{
  %% The code I acutally used to get the averages

  news = textread(['good_' num2str(num) '.txt'], '%s');

  imgs = cell(length(news), 1);
  signals = cell(length(news), 1);
  times = NaN(length(news), 1);

  for i=1:length(news)
    load(news{i});

    [fraction, max_width, cell_width, img] = domain_expansion(mymovie, opts);
    pos = [1:size(img,1)].';
    center = (size(img,2)-1)/2+1;
    indx = find(isnan(fraction), 1, 'last');
    if (isempty(indx))
      indx = 1;
    else
      fraction(indx) = 0;
    end
    poly = [[center-2*max_width*fraction(indx:end) pos(indx:end)]; [center+2*max_width*fraction(end:-1:indx) pos(end:-1:indx)]];
    mask = roipoly(size(img,1),size(img,2), poly(:,1), poly(:,2));
    weight = exp(-bwdist(mask)/(2*5^2));

    img2 = img.*weight;

    t = get_manual_timing(mymovie, opts);

    signals{i} = img;
    imgs{i} = img2;
    times(i) = t(1);
  end

  uuid = num2str(now+cputime);

  signals_full = cell(10,3);
  for i=1:10
    [signals_full{i,1}, signals_full{i,2}, signals_full{i,3}] = simultaneous_registration(imgs, times);
    save(['aligned_domains_' num2str(num) '-' uuid '.mat'], 'signals_full', 'imgs', 'signals');
  end
  %}

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
