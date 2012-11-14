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
    domain = domain(1:time(end), :);

    size_diff = size(signals, 2) - size(domain, 2);
    if (size_diff > 0)
      domain = padarray(domain, [0 size_diff/2], NaN);
    end
    if (~isempty(signals) & size_diff < 0)
      signals = padarray(signals, [0 -size_diff/2 0], NaN);
    end

    switch (sync_type)
      case 'fraction'
        align_time = find(fraction >  sync_param, 1, 'first');
        if (isempty(align_time))
          align_time = time(end);
        end
      case 'lsr'
        domain = domain / mymean(domain(:));
        [signals, align_time] = find_min_residue(signals, align_time, domain, time(1));
    end
  end

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

function [signal, indx] = find_min_residue(signal, indx, candidate, cand_indx)

  if (isempty(signal))
    signal = candidate;
    indx = cand_indx;

    return;
  end
  
  [h,w,n] = size(signal);
  npts = size(candidate, 1);
  cand_pos = 1:npts;
  sign_pos = [1:h];
  rel_pos = sign_pos - indx;

  fval = mymean(mymean(candidate(1:5,:), 1), 2);
  lval = mymean(mymean(candidate(end-5:end, :), 1), 2);

  syncs = [-15:15];
  nsync = length(syncs);
  residues = NaN(nsync, 1);

  for i=1:nsync
    tmp_pos = cand_pos - syncs(i) - cand_indx;
    goods = (tmp_pos >= rel_pos(1) & tmp_pos <= rel_pos(end));
    sign_goods = (rel_pos >= tmp_pos(1) & rel_pos <= tmp_pos(end));
    
    if (sum(goods)/npts < 0.5 | sum(sign_goods)/h < 0.5)
      residues(i) = Inf;
    else
      tmp_vals = NaN(h,w);
      tmp_vals(sign_pos(sign_goods), :) = candidate(cand_pos(goods), :);
      tmp_vals(rel_pos < tmp_pos(1), :) = fval;
      tmp_vals(rel_pos > tmp_pos(end), :) = lval;
      %ends = (rel_pos > tmp_pos(end));
      %tmp_vals(ends, :) = repmat(lval, sum(ends), 1);

      [means, stds] = mymean(cat(3, signal, tmp_vals), 3);
      residues(i) = mymean(stds(:));
    end
  end

  [val, rel_indx] = min(residues);

  tmp_pos = cand_pos - syncs(rel_indx) - cand_indx;
  ends = (tmp_pos > rel_pos(end));

  if (any(ends))
    nends = sum(ends);
    signal = padarray(signal, [nends 0 0], NaN, 'post');

    [h,w,n] = size(signal);
    sign_pos = [1:h];
    rel_pos = sign_pos - indx;
  end

  goods = (tmp_pos >= rel_pos(1) & tmp_pos <= rel_pos(end));
  sign_goods = (rel_pos >= tmp_pos(1) & rel_pos <= tmp_pos(end));
  
  tmp_vals = NaN(h,w);
  tmp_vals(sign_pos(sign_goods), :) = candidate(cand_pos(goods), :);

  signal = cat(3, signal, tmp_vals);

  return;
end
