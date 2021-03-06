function [embryos] = split_cells(embryos_estim, angle_thresh, max_ratio, max_distance, max_score, max_overlap)

  if (~iscell(embryos_estim))
    embryos_estim = {embryos_estim};
  end

  nframes = length(embryos_estim);

  embryos = cell(nframes, 1);

  for nimg = 1:nframes

    estim = embryos_estim{nimg};
    if (isempty(estim))
      continue;
    end

    pts = find(isnan(estim(:,1)));
    if (isempty(pts))
      pts = size(estim, 1);
    end
    pts = [1; pts];

    for i=1:length(pts)-1
      tmp_estim = estim(pts(i):pts(i+1), :);
      tmp_estim = tmp_estim(~any(isnan(tmp_estim), 2), :);

      [pac, indxs] = impac(tmp_estim);

      concaves = compute_concavity(pac, angle_thresh);

      ellipses = fit_ellipse_segments(tmp_estim, indxs(concaves), max_ratio, max_distance, max_score, max_overlap);

      if (i == 1)
        embryos{nimg} = ellipses;
      else
        embryos{nimg} = [embryos{nimg}; ellipses];
      end
    end
  end

  if (nframes == 1)
    embryos = embryos{1};
  end

  return;
end
