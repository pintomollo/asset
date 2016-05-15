function embryos = filter_embryos(embryos, angle_thresh, max_ratio, max_distance, max_score, max_overlap, max_area_diff)

  if (nframes == 1)
    all_ellipses = all_ellipses{1};
  elseif (max_nellipses > 1)
    display('There seems to be several cells in the images, double checking if everything is fine !');

    dist_thresh = (max_distance / 3)^2;
    start_index = nframes;
    for i=1:nframes
      if (~isempty(all_ellipses{i}))
        real_ell = NaN(size(all_ellipses{i},1), ndata, nframes);
        real_ell(:,:,1) = all_ellipses{i};
        start_index = i;

        break;
      end
    end

    for i=start_index+1:nframes
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
          real_ell(end+1,:,:) = NaN;
          real_ell(end, :, 1) = tmp_ell(k,:);
        end
      end
    end

    goods = sum(~isnan(real_ell),3);
    goods = (goods(:,1) > nframes/3);

    avg_ell = NaN(0,ndata);

    for i=1:size(real_ell, 1)
      if (goods(i))
        tmp = real_ell(i, :, :);
        tmp = reshape(tmp(~isnan(tmp)), 5, []);
        avg_ell(end+1,:) = median(tmp, 2);

        %draw_ellipse(avg_ell(end, 1:2), avg_ell(end, 3:4), avg_ell(end, 5), 'r');
      end
    end
    avg_area = prod(avg_ell(:,3:4), 2)*pi;
    mean_area = mean(avg_area);
    goods = (avg_area <= avg_area * max_area_diff & ...
             avg_area >= avg_area / max_area_diff);

    avg_ell = avg_ell(goods, :);
    ntargets = size(avg_ell, 1);

    for nimg = 1:nframes
      
      estim = all_estim{nimg};
      
      if (isempty(estim))
        continue;
      end
          
      [pac, indxs] = impac(estim);
      concaves = compute_concavity(pac, angle_thresh);

      sides = (any(estim == 2 | bsxfun(@eq, estim, imgsize-1), 2));
      borders = (xor(sides, sides([2:end 1])));
      sides(borders) = false;

      borders = borders | isnan(estim(:,1));
      borders(indxs(concaves)) = true;
      
      labels = cumsum(double(borders), 1);

      estim = estim(~sides, :);
      labels = labels(~sides);

      nlabel = hist(labels, labels(end));
      label_indx = [1:length(nlabel)];

      too_few = ismember(labels, label_indx(nlabel < 10));
      estim = estim(~too_few, :);
      labels = labels(~too_few);

      indxs = unique(labels);

      ellipses = NaN(ntargets, 5);
      scores = Inf(ntargets, 1);

      for i = 1:ntargets
        %current = ismember(labels, indxs);
        %tmp_pts = estim(current, :);
        ell_pts = carth2elliptic(estim, avg_ell(i, 1:2), avg_ell(i, 3:4), avg_ell(i, 5), 'radial');
        dist = abs(ell_pts(:,2) - 1);

        valids = (dist <= 2*max_score);
        npts = mymean(double(valids), 1, labels);

        fit_indxs = indxs(npts > 0.333);
        current = ismember(labels, fit_indxs);
        [ellipses(i, :)] = fit_distance(estim(current, :));
        scores(i) = overlaps(avg_ell(i,:), ellipses(i,:), NaN); 
      end

      goods = (scores > 0.75 & (ellipses(:, 4) ./ ellipses(:, 3)) >= max_ratio);
      ellipses(~goods, :) = avg_ell(~goods, :);

       % for j = 1:size(ellipses, 1)
       %   draw_ellipse(ellipses(j, 1:2), ellipses(j, 3:4), ellipses(j, 5));
       % end

        all_ellipses{nimg} = ellipses;
    end
  end

  return;
end
