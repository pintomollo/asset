function [signal, indx, rel_indx] = find_min_residue(signal, indx, candidate, cand_indx)

  if (isempty(signal))
    signal = candidate;
    indx = cand_indx;
    rel_indx = 0;

    return;
  end
  
  [h,w,n] = size(signal);
  npts = size(candidate, 1);
  cand_pos = 1:npts;
  sign_pos = [1:h];
  rel_pos = sign_pos - indx;

  fval = mymean(mymean(candidate(1:5,:), 1), 2);
  lval = mymean(mymean(candidate(end-5:end, :), 1), 2);

  syncs = ceil(size(candidate, 1) / 4);
  syncs = [-syncs:syncs];
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
  rel_indx = syncs(rel_indx) + cand_indx;

  return;
end
