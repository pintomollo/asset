function embryos = filter_embryos(all_ellipses, all_estims, frame_indx, slice_indx, angle_thresh, max_ratio, max_distance, max_score, max_overlap, max_area_diff)

  nimgs = length(all_ellipses);
  nframes = max(frame_indx);
  nslices = max(slice_indx);

  for i=1:nimgs
    if (~isempty(all_ellipses{i}))
      all_ellipses{i}(:,end+1) = frame_indx(i);
      all_ellipses{i}(:,end+1) = slice_indx(i);
    end
  end

  nembryos = cellfun(@(x)size(x,1), all_ellipses);
  all_embryos = cat(1, all_ellipses{:});

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
    indxs = new_indxs(indxs);
  end

  frames = all_embryos(:,end-1);
  slices = all_embryos(:,end);

  ntargets = max(nembryos);
  nframes = max(frames);
  nslices = max(slices);

  full_indxs = (frames-1)*ntargets*nslices + (slices-1)*ntargets + clusts;

  new_embryos = NaN(nframes*nslices*ntargets, 5);

  new_embryos(full_indxs, :) = all_embryos(:, 1:5);
  new_embryos = permute(reshape(new_embryos.', [5 ntargets nslices nframes]), [2 1 3 4]);

  avg_ell = mymean(new_embryos, 4);
  goods = (sum(~any(isnan(new_embryos), 2), 4) > nframes/3);

  keyboard



  %{
  avg_ell = mymean(all_embryos, 1, clusts);
  nembryos = length(indxs);

  avg_area = prod(avg_ell(:,3:4), 2)*pi;
  mean_area = mean(avg_area);
  goods = (avg_area <= avg_area * max_area_diff & ...
           avg_area >= avg_area / max_area_diff);

  avg_ell = avg_ell(goods, :);
  ntargets = size(avg_ell, 1);
  real_ell = NaN(ntargets, 5, nframes);
  %}

  %full_indxs = (frames-1)*nembryos + clusts;
  %new_embryos = NaN(nframes*nembryos, 5);

  %new_embryos(full_indxs, :) = all_embryos(:, 1:5);
  %new_embryos2 = permute(reshape(new_embryos.', [5 nembryos nframes]), [2 1 3]);

  %keyboard
  %%%

  %{
  last_indx = 0;
  for i=1:nframes
    tmp_ell = all_ellipses{i};
    for k=1:size(tmp_ell, 1)
      found = false;
      mpos = mymean(real_ell(:,1:2,:), 3);
      for j=1:size(mpos, 1)
        if (sum((mpos(j, :) - tmp_ell(k, 1:2)).^2) < dist_thresh)
          found = true;
          real_ell(j, :, i) = tmp_ell(k, :);
          break;
        end
      end

      if (~found)
        last_indx = last_indx + 1;
        real_ell(last_indx, :, 1) = tmp_ell(k,:);
      end
    end
  end

  goods = sum(~isnan(real_ell),3);
  goods = (goods(:,1) > nframes/3);
  real_ell = real_ell(goods,:,:);

  avg_ell = median(real_ell, 3, 'omitnan');
%  for i=1:size(real_ell, 1)
%    tmp = real_ell(i, :, :);
%    tmp = reshape(tmp(~isnan(tmp)), 5, []);
%    avg_ell(end+1,:) = median(tmp, 2);
%  end
  avg_area = prod(avg_ell(:,3:4), 2)*pi;
  mean_area = mean(avg_area);
  goods = (avg_area <= avg_area * max_area_diff & ...
           avg_area >= avg_area / max_area_diff);

  avg_ell = avg_ell(goods, :);
  ntargets = size(avg_ell, 1);
  %}

  for nimg = 1:nframes

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

    ellipses = NaN(ntargets, 5);
    scores = Inf(ntargets, 1);

    for i = 1:ntargets
      ell_pts = carth2elliptic(estim, avg_ell(i, 1:2), avg_ell(i, 3:4), avg_ell(i, 5), 'radial');
      dist = abs(ell_pts(:,2) - 1);

      valids = (dist <= 2*max_score);
      npts = mymean(double(valids), 1, labels);

      fit_indxs = indxs(npts > 0.333);
      current = ismember(labels, fit_indxs);
      current_concaves = ismember(concave_labels, fit_indxs);
      new_ellipse = fit_ellipse_segments(estim(current,:), concave_labels(current_concaves), max_ratio, max_distance, max_score, max_overlap);

      if (isempty(new_ellipse))
        real_ell(i,:,nimg) = NaN;
      else
        real_ell(i,:,nimg) = new_ellipse;
      end
    end
  end

  %% NEED TO SETUP THE STUFF IN SHAPE & INTERPOLATE FOR NAN ?
  embryos = real_ell;

  return;
end
