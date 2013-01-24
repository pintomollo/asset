function [signal, indx, rel_indx] = find_min_residue(signal, indx, candidate, cand_indx, thresh, syncs)

  if (isempty(signal))
    signal = candidate;
    indx = cand_indx;
    rel_indx = 0;

    return;
  end

  if (nargin < 5)
    thresh = 0.5;
    syncs = -1;
  elseif (nargin < 6)
    if (thresh > 1)
      syncs = thresh;
      thresh = 0.5;
    else
      syncs = -1;
    end
  end
  
  [w,h,n] = size(signal);
  [npos, npts, nlayers] = size(candidate);

  sstd = std(signal(~isnan(signal)));
  cstd = std(candidate(~isnan(candidate)));

  if (w ~= npos)
    [X, Y] = meshgrid([1:npts], 1+([0:npos-1]*(w-1)/(npos-1)).');
    candidate = bilinear_mex(candidate, X, Y, [2 2]);
    [npos, npts] = size(candidate);
  end

  cand_pos = 1:npts;
  sign_pos = [1:h];
  rel_pos = sign_pos - indx;

  fval = mymean(mymean(candidate(1:5,:), 2), 1);
  lval = mymean(mymean(candidate(end-5:end, :), 2), 1);

  if (syncs < 0)
    syncs = ceil(size(candidate, 2) / 4);
  end

  syncs = [-syncs:syncs].';
  nsync = length(syncs);
  residues = NaN(nsync, 3);

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
      residues(i,1) = mymean(stds(:));

%      if (all(goods))
%        residues(i,2)  = cstd;
%      else
        tmp_val = candidate(:, cand_pos(~goods));
        [junk, residues(i,2)] = mymean(tmp_val(:));
%      end
%      if (all(sign_goods))
%        residues(i,3)  = sstd;
%      else
        tmp_val = signal(:,sign_pos(~sign_goods),:);
        [junk, residues(i,3)] = mymean(tmp_val(:));
%      end
    end
  end

  bads = isnan(residues(:,2));
  goods = ~bads;
  if (any(goods))
    findx = find(goods, 1, 'first');
    lindx = find(goods, 1, 'last');
    bads([1:findx lindx:end]) = false;
    residues(bads, 2) = interp1q(syncs(goods), residues(goods, 2), syncs(bads));
  end

  bads = isnan(residues(:,3));
  goods = ~bads;
  if (any(goods))
    findx = find(goods, 1, 'first');
    lindx = find(goods, 1, 'last');
    bads([1:findx lindx:end]) = false;
    residues(bads, 3) = interp1q(syncs(goods), residues(goods, 3), syncs(bads));
  end

  if (none)
    rel_indx = 0;
    tmp_pos = cand_pos - cand_indx;
    ends = (tmp_pos > rel_pos(end));
    starts = (tmp_pos < rel_pos(1));

    indx = 0;

    warning('Signal could not be aligned !');
  else
    [val, rel_indx] = min(residues(:,1));
    tmp_pos = cand_pos - syncs(rel_indx) - cand_indx;
    starts = (tmp_pos < rel_pos(1));
    ends = (tmp_pos > rel_pos(end));

    rel_indx = syncs(rel_indx) + cand_indx;
  end

  nstarts = sum(starts);
  if (nstarts > 0)
    signal = padarray(signal, [0 nstarts 0], NaN, 'pre');
  end

  if (any(ends))
    nends = sum(ends);
    signal = padarray(signal, [0 nends 0], NaN, 'post');
  end

  [w,h,n] = size(signal);
  sign_pos = [1:h];
  rel_pos = sign_pos - indx - nstarts;

  goods = (tmp_pos >= rel_pos(1) & tmp_pos <= rel_pos(end));
  sign_goods = (rel_pos >= tmp_pos(1) & rel_pos <= tmp_pos(end));
  
  tmp_vals = NaN(w,h);
  tmp_vals(:, sign_pos(sign_goods)) = candidate(:, cand_pos(goods));
%  tmp_vals(:, sign_pos(sign_goods)) = candidate;

  signal = cat(3, signal, tmp_vals);

  return;
end
