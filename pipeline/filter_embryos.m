function embryos = filter_embryos(all_ellipses, all_estims, angle_thresh, max_ratio, max_distance, max_score, max_overlap, max_area_diff)

  nframes = length(all_ellipses);
  nembryos = cellfun(@(x)size(x,1), all_ellipses);
  real_ell = NaN(max(nembryos), 5, nframes);

  last_indx = 0;
  dist_thresh = (max_distance / 3)^2;
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
      %[ellipses(i, :)] = fit_distance(estim(current, :));
      %scores(i) = overlaps(avg_ell(i,:), ellipses(i,:), NaN); 
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
