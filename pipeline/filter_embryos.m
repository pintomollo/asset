function embryos = filter_embryos(all_ellipses, all_estims, frame_indx, slice_indx, angle_thresh, max_ratio, max_distance, max_score, max_overlap, max_area_diff)

  nimgs = length(all_ellipses);
  nframes = max(frame_indx);

  for i=1:nimgs
    if (~isempty(all_ellipses{i}))
      all_ellipses{i}(:,end+1) = frame_indx(i);
      all_ellipses{i}(:,end+1) = slice_indx(i);
    end
  end

  nembryos = cellfun(@(x)size(x,1), all_ellipses);
  all_embryos = cat(1, all_ellipses{:});
  nparams = size(all_embryos, 2) - 2;

  dist_thresh = (max_distance / 3);
  clusts = cluster_distance(all_embryos(:,1:2), dist_thresh);

  indxs = [1:max(nembryos)];
  counts = hist(clusts, indxs);
  indxs = indxs(counts > nframes/3);
  goods = ismember(clusts, indxs);

  all_embryos = all_embryos(goods, :);

  clusts = clusts(goods);
  if (max(indxs)~=length(indxs))
    new_indxs = [1:max(indxs)];
    new_indxs(indxs) = [1:length(indxs)];

    clusts = new_indxs(clusts);
  end

  frames = all_embryos(:,end-1);
  slices = all_embryos(:,end);

  ntargets = max(nembryos);
  nframes = max(frames);
  nslices = max(slices);

  full_indxs = (frames-1)*ntargets*nslices + (slices-1)*ntargets + clusts;

  embryos = NaN(nframes*nslices*ntargets, nparams);

  embryos(full_indxs, :) = all_embryos(:, 1:nparams);
  embryos = permute(reshape(embryos.', [nparams ntargets nslices nframes]), [2 1 3 4]);

  goods = (sum(~any(isnan(embryos), 2), 4) > nframes/3);
  all_goods = (sum(goods, 3) > nslices/2);

  embryos = embryos(all_goods,:,:,:);
  goods = goods(all_goods,:,:);

  if (~any(goods(:)))
    embryos = [];
    return;
  end

  ntargets = size(embryos, 1);

  avg_ell = median(embryos, 4, 'omitnan');
  avg_area = prod(avg_ell(:,3:4,:), 2)*pi;
  mean_area = mean(avg_area(:));

  area_goods = (avg_area <= mean_area * max_area_diff & ...
           avg_area >= mean_area / max_area_diff);
  goods = goods & area_goods;

  embryos(repmat(~goods, [1 nparams 1 nframes])) = NaN;

  avg_ell = median(embryos, 4, 'omitnan');
  avg_avg = mymean(avg_ell, 3);

  if (any(isnan(avg_avg(:))))
    embryos = [];
    return;
  end

  for nimg = 1:nimgs

    estim = all_estims{nimg};

    if (isempty(estim))
      continue;
    end

    [pac, indxs] = impac(estim);
    concaves = compute_concavity(pac, angle_thresh);

    borders = isnan(estim(:,1));
    borders(indxs(concaves)) = true;

    labels = cumsum(double(borders), 1);
    concave_labels = labels(indxs(concaves));

    nlabel = hist(labels, labels(end));
    label_indx = [1:length(nlabel)];

    too_few = ismember(labels, label_indx(nlabel < 10));
    estim = estim(~too_few, :);
    labels = labels(~too_few);

    indxs = unique(labels);

    for i = 1:ntargets

      curr_avg = avg_ell(i, :, slice_indx(nimg));
      if (isnan(curr_avg))
        curr_avg = avg_avg(i, :);
      end

      ell_pts = carth2elliptic(estim, curr_avg(1:2), curr_avg(3:4), curr_avg(5), 'radial');
      dist = abs(ell_pts(:,2) - 1);

      valids = (dist <= 2*max_score);
      npts = mymean(double(valids), 1, labels);

      fit_indxs = indxs(npts > 0.333);
      current = ismember(labels, fit_indxs);
      current_concaves = ismember(concave_labels, fit_indxs);
      ellipse = fit_ellipse_segments(estim(current,:), ...
                                     concave_labels(current_concaves), max_ratio, ...
                                     max_distance, max_score, max_overlap);

      if (isempty(ellipse))
        embryos(i,:,slice_indx(nimg), frame_indx(nimg)) = NaN;
      else
        embryos(i,:,slice_indx(nimg), frame_indx(nimg)) = ellipse;
      end
    end
  end

  return;
end
