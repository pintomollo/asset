function [all_pts, pos] = domain_scaling(fnames, signal_type, sync_type)
  
  tmp_name = strrep(fnames, '*', 'domains');

  if (nargin == 1)
    signal_type = 'width';
    sync_type = 'fraction';
  end

  if (exist(tmp_name, 'file'))
    load(tmp_name);
  else

    movies = dir(fnames);
    nmovies = length(movies);

    datas = cell(nmovies, 1);
    widths = NaN(nmovies, 1);
    lengths = NaN(nmovies, 1);

    for i=1:nmovies
      load(movies(i).name);

      opts.recompute = false;
      [fraction, width, cell_length] = domain_expansion(mymovie, opts);

      datas{i} = fraction;
      widths(i) = width;
      lengths(i) = cell_length;
    end

    save(tmp_name, 'datas', 'widths', 'lengths');
  end

  nmovies = length(datas);
  all_pts = NaN(0,3);

  shuffle = randperm(nmovies);

  for i=shuffle
    fraction = datas{i};
    fraction(isnan(fraction)) = 0;

    switch (sync_type)
      case 'fraction'
        pt = find(fraction > 0.775, 1);
    end

    switch (signal_type)
      case 'width'
        fraction = fraction * widths(i);
      case 'cell'
        fraction = fraction * 2*widths(i)./lengths(i);
    end

    indexes = [1:length(fraction)].';
    candidate = [fraction, ones(size(fraction))*i, (indexes)];

    switch (sync_type)
      case 'lsr'
        pt = find_min_residue(all_pts, candidate);
    end

    candidate(:,end) = candidate(:,end) - pt;

    all_pts = extend_signal(all_pts, candidate);
  end

  all_pts(isnan(all_pts(:,1))) = 0;
  all_pts(:,end)  = all_pts(:,end)*10;

  all_pts = all_pts(:, [2 3 1]);
  all_pts = sortrows(all_pts);

  pos = unique(all_pts(:,2));
  all_pts = reshape(all_pts(:,3), length(pos), []);

  %[avg, stds] = mymean(all_pts(:,1).*widths(all_pts(:,2)), 1, all_pts(:,3));

  %figure;scatter(all_pts(:,3), all_pts(:,1).*widths(all_pts(:,2)), 'b');
  %hold on;
  %plot(pts, avg, 'k');
  %plot(pts, avg+stds, 'r')
  %plot(pts, avg-stds, 'r')

  %{
  [avg, stds] = mymean(all_pts(:,1), 1, all_pts(:,3));

  figure;scatter(all_pts(:,3), all_pts(:,1), 'b');
  hold on;
  plot(pts, avg, 'k');
  plot(pts, avg+stds, 'r')
  plot(pts, avg-stds, 'r')

  title([tmp_name ':' signal_type '(' sync_type ')']);

  all_speeds = NaN(0, 2);
  for i=1:nmovies

    signal = all_pts(all_pts(:, 2) == i, :);
    signal = sortrows(signal(:,[3 1 2]));
    speed = differentiator(signal(:,1), signal(:,2), 'replicate');

    all_speeds = [all_speeds; [speed signal(:, 1)]];
  end

  [avg, stds] = mymean(all_speeds(:,1), 1, all_speeds(:,2));

  figure;scatter(all_speeds(:,2), all_speeds(:,1), 'b');
  hold on;
  plot(pts, avg, 'k');
  plot(pts, avg+stds, 'r')
  plot(pts, avg-stds, 'r')

  title([tmp_name ':' signal_type '(' sync_type ')']);

  %}

%  ratios = 2*widths./lengths;

%  [avg, stds] = mymean(all_pts(:,1).*ratios(all_pts(:,2)), 1, all_pts(:,3));

%  figure;scatter(all_pts(:,3), all_pts(:,1).*ratios(all_pts(:,2)), 'b');
%  hold on;
%  plot(pts, avg, 'k');
%  plot(pts, avg+stds, 'r')
%  plot(pts, avg-stds, 'r')

  return;
end

function signal = extend_signal(signal, candidate)

  if (isempty(signal))
    signal = candidate;

    return;
  end

  pos = unique(signal(:,end));
  cand_pos = unique(candidate(:, end));
  signals = unique(signal(:,end-1));

  ends = [cand_pos(end)+1:pos(end)].';
  starts = [pos(1):cand_pos(1)-1].';
  candidate = [[zeros(size(starts)), candidate(1,2)*ones(size(starts)), starts]; ...
              candidate; ...
              [candidate(end,1)*ones(size(ends)), candidate(1,2)*ones(size(ends)), ends]];

  firsts = signal(signal(:,end)==pos(1),:);
  firsts(:,1) = 0;
  for i=cand_pos(1):pos(1)-1
    firsts(:,end) = i;
    signal = [signal; firsts];
  end

  lasts = signal(signal(:,end)==pos(end),:);
  for i=pos(end)+1:cand_pos(end)
    lasts(:,end) = i;
    signal = [signal; lasts];
  end

  signal = [signal; candidate];

  return;
end

function indx = find_min_residue(signal, candidate)

  if (isempty(signal))
    indx = ceil(size(candidate, 1) / 2);

    return;
  end
  
  pos = unique(signal(:,end));
  npts = size(candidate, 1);

  cand_pos = unique(candidate(:, end));
  fval = 0;
  lval = mean(candidate(end-5:end, 1));

  residues = NaN(npts, 1);

  for i=1:npts
    tmp_pos = cand_pos - i;
    tmp_vals = NaN(size(pos));
    tmp_vals(pos < tmp_pos(1)) = fval;
    tmp_vals(pos > tmp_pos(end)) = lval;

    [means, stds] = mymean([signal(:,1); candidate(:,1); tmp_vals], 1, [signal(:,end); tmp_pos; pos]);
    residues(i) = mymean(stds);
  end

  [val, indx] = min(residues);

  return;
end
