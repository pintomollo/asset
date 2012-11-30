function [signal, indx, rel_indx] = find_min_residue(signal, indx, candidate, cand_indx, thresh)

  if (isempty(signal))
    signal = candidate;
    indx = cand_indx;
    rel_indx = 0;

    return;
  end

  if (nargin < 5)
    thresh = 0.5;
  end
  
  [w,h,n] = size(signal);
  npts = size(candidate, 2);
  cand_pos = 1:npts;
  sign_pos = [1:h];
  rel_pos = sign_pos - indx;

  fval = mymean(mymean(candidate(1:5,:), 2), 1);
  lval = mymean(mymean(candidate(end-5:end, :), 2), 1);

  syncs = ceil(size(candidate, 2) / 4);
  syncs = [-syncs:syncs];
  nsync = length(syncs);
  residues = NaN(nsync, 1);

  cand_thresh = min(thresh, npts/h);
  sign_thresh = min(thresh, h/npts);

  none = true;
  for i=1:nsync
    tmp_pos = cand_pos - syncs(i) - cand_indx;
    goods = (tmp_pos >= rel_pos(1) & tmp_pos <= rel_pos(end));
    sign_goods = (rel_pos >= tmp_pos(1) & rel_pos <= tmp_pos(end));
    
    if (sum(goods)/npts < sign_thresh | sum(sign_goods)/h < cand_thresh)
      residues(i) = Inf;
    else
      none = false;

      tmp_vals = NaN(w,h);
      tmp_vals(:, sign_pos(sign_goods)) = candidate(:, cand_pos(goods));
      tmp_vals(:, rel_pos < tmp_pos(1)) = fval;
      tmp_vals(:, rel_pos > tmp_pos(end)) = lval;
      %ends = (rel_pos > tmp_pos(end));
      %tmp_vals(ends, :) = repmat(lval, sum(ends), 1);

      [means, stds] = mymean(cat(3, signal, tmp_vals), 3);
      residues(i) = mymean(stds(:));
    end
  end

  if (none)
    rel_indx = 0;
    tmp_pos = cand_pos - cand_indx;
    ends = (tmp_pos > rel_pos(end));

    warning('Signal could not be aligned !');
  else
    [val, rel_indx] = min(residues);
    tmp_pos = cand_pos - syncs(rel_indx) - cand_indx;
    ends = (tmp_pos > rel_pos(end));

    rel_indx = syncs(rel_indx) + cand_indx;
  end

  if (any(ends))
    nends = sum(ends);
    signal = padarray(signal, [0 nends 0], NaN, 'post');

    [w,h,n] = size(signal);
    sign_pos = [1:h];
    rel_pos = sign_pos - indx;
  end

  goods = (tmp_pos >= rel_pos(1) & tmp_pos <= rel_pos(end));
  sign_goods = (rel_pos >= tmp_pos(1) & rel_pos <= tmp_pos(end));
  
  tmp_vals = NaN(w,h);
  tmp_vals(:, sign_pos(sign_goods)) = candidate(:, cand_pos(goods));

  signal = cat(3, signal, tmp_vals);

  return;
end
