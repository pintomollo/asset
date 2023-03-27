function display_angular(field, correction, time_corr)

  load('results_segment');

  if (nargin < 2 & strncmp(mymovie.experiment, 'GZ920_MERGE_RECOS-', 18))
    correction = false(1,75);
    correction([1:15 46:60]) = true;
    time_corr = ones(1,75);
    time_corr(1:15) = 690.162;
    time_corr(16:30) = 1040.16;
    time_corr(31:45) = 1120.163;
    time_corr(46:60) = 950.158;
    time_corr(61:75) = 1010.174;
  end
  if (nargin == 0)
    field = 'dic';
  end

  %orig_errors = mymovie2.(field).errors;
  %estim_errors = orig_errors(:,:,1,:);
  %orig_errors = orig_errors(:,:,2,:);
  optim_errors = mymovie.(field).errors;
  track_errors = trackings.(field).errors;
  %intra_errors = [];
  %for i=1:length(trackings.(field).child)
  %  intra_errors = cat(3, intra_errors, trackings.(field).child(i).errors);
  %end

  times = frame_timing(mymovie.eggshell);
  times = times - time_corr;

  [b,indx] = histc(times, [-1020:60:180]);

  %orig_errors(:,correction,:,:) = orig_errors(:,correction,:,[(end/2)+1:end 1:end/2],:);
  %estim_errors(:,correction,:,:) = estim_errors(:,correction,:,[(end/2)+1:end 1:end/2],:);
  optim_errors(:,correction,:,:) = optim_errors(:,correction,:,[(end/2)+1:end 1:end/2],:);
  %track_errors(:,correction,:,:) = track_errors(:,correction,:,[(end/2)+1:end 1:end/2],:);
  %intra_errors(:,correction,:,:) = intra_errors(:,correction,:,[(end/2)+1:end 1:end/2],:);
  
  %diff_errors = repmat(optim_errors,[1 1 size(track_errors,3) 1]) - track_errors;

  [egg, cortex, index_egg, index_cortex] = compute_errors(optim_errors, 4);
  %[egg, cortex, index_egg, index_cortex] = compute_errors(diff_errors, 4);
  all = [egg;cortex];
  all_index = [index_egg; index_cortex+max(index_egg)];

  [errs, stds] = mymean(all, [], all_index);
  errs = [errs(1:end/2) errs(end/2+1:end)];
  stds = [stds(1:end/2) stds(end/2+1:end)];
  
  [H,p] = myttest([egg;egg],[ones(size(index_egg)); index_egg+1]);
  [H(:,1) p(:,1)].'
  [H,p] = myttest([cortex;cortex],[ones(size(index_cortex)); index_cortex+1]);
  [H(:,1) p(:,1)].'

  figure;
  subplot(211);
  barweb(errs * 100, stds * 100, [], '', '', '', 'Angular',[],'y',[],1);

  %[egg, cortex, index_egg, index_cortex] = compute_errors(diff_errors, 2);
  [egg, cortex, index_egg, index_cortex] = compute_errors(optim_errors, 2);

  for i=1:length(indx)
    index_egg(index_egg==i) = indx(i);
    index_cortex(index_cortex==i) = indx(i);
  end

  all = [egg;cortex];
  all_index = [index_egg; index_cortex+max(index_egg)];
  groups_index = [ones(size(index_egg)); 2*ones(size(index_cortex))];

  [H,p] = myttest([egg;egg],[ones(size(index_egg)); index_egg+1]);
  [H(:,1) p(:,1)].'
  [H,p] = myttest([cortex;cortex],[ones(size(index_cortex)); index_cortex+1]);
  [H(:,1) p(:,1)].'

  [errs, stds] = mymean(all, [], all_index);
  errs = [errs(1:end/2) errs(end/2+1:end)];
  stds = [stds(1:end/2) stds(end/2+1:end)];
  
  %errs = reshape(errs, [], 2);
  %stds = reshape(stds, [], 2);
  subplot(212);
  barweb(errs * 100, stds * 100, [], '', '', '', 'Time',[],'y',[],1);

  return;
end

function [eggs, cors, indxe, indxc] = compute_errors(errors, dim, eggs, cors, indxe, indxc)

  if (nargin == 2)
    eggs = [];
    cors = [];
    indxe = [];
    indxc = [];
    shift = 0;
  else
    shift = max(indxe);
  end

  for i=1:size(errors,dim)
    switch dim
      case 2
        egg = errors(1,i,:,:);
        cor = errors(2,i,:,:);
      case 3
        egg = errors(1,:,i,:);
        cor = errors(2,:,i,:);
      case 4
        egg = errors(1,:,:,i);
        cor = errors(2,:,:,i);
    end

    egg = egg(:);

    cor = cor(:);

    eggs = [eggs; egg];
    cors = [cors; cor];

    indxe = [indxe; ones(size(egg))*(shift + i)];
    indxc = [indxc; ones(size(cor))*(shift + i)];
  end

  return;
end
