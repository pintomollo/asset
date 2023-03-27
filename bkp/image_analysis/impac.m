function [pac, hits] = impac(contour, nsteps, thresh, method)
  % Wang, W. X. (1998). Binary Image Segmentation of Aggregates Based on Polygonal Approximation and Classification of Concavities. Pattern Recognition, 31(10), 1503–1524. doi:10.1016/S0031-3203(97)00145-3

  switch nargin
    case 1
      nsteps = [];
      thresh = [];
      method = 'recursive';
    case 2
      if (ischar(nsteps))
        method = nsteps;
        thresh = [];
      else
        thresh = nsteps;
        method = 'recursive';
      end
      nsteps = [];
    case 3
      if (ischar(thresh))
        method = thresh;

        switch method
          case 'recursive'
            thresh = nsteps;
            nsteps = [];
          otherwise
            thresh = [];
        end
      else
        method = 'iterative';
      end
  end

  is_repeated = true;
  if (all(contour(1,:) ~= contour(end, :)))
    contour = contour([1:end 1], :);
    is_repeated = false;
  end
  npts = size(contour, 1);

  if (isempty(thresh))
    % Rough incircle radius estimation (from the bounding box)
    thresh = 0.1 * min((max(contour, [], 1) - min(contour, [], 1))) / 2.5;

    if (method(1) == 'i')
      thresh = thresh / 4;
    end
  end

  if (isempty(nsteps))
    nsteps = ceil(npts / 10);
  end
  nsteps = nsteps + 1;

  pac = NaN(size(contour));
  switch method
    case 'recursive'
      if (hypot(contour(end-1, 1) - contour(end, 1), contour(end-1,2)-contour(end,2)) > thresh)
        indx = size(contour, 1)-1;
      else
        indx = ceil(size(contour, 1) / 2);
      end

      if (indx > 0)
        pac(1:indx, :) = pac_recursive(contour(1:indx, :), thresh);
      end
      pac(indx:end, :) = pac_recursive(contour(indx:end, :), thresh);
    case 'iterative'
      pac = pac_iterative(contour, nsteps, thresh);
  end

  if (~is_repeated && all(pac(1,:) == pac(end, :)))
    pac = pac([1:end-1], :);
  end
  hits = find(~any(isnan(pac), 2));
  pac = pac(hits, :);

  %figure;
  %scatter(contour(:,1), contour(:,2))
  %hold on;
  %scatter(pac(:,1), pac(:,2), 'r');

  return;
end

function pac = pac_recursive(contour, thresh)

  pac = NaN(size(contour));

  end_indx = size(pac, 1) + 1;
  segment = zeros(1, 2);

  while (all(segment == 0) & end_indx > 1)
    end_indx = end_indx - 1;
    segment = contour(end_indx, :) - contour(1, :);
  end

  pac([1 end_indx], :) = contour([1 end_indx], :);

  if (end_indx < 3 | size(contour, 1) < 3)
    return;
  end

  dist = abs(segment(1, 1)*(pac(1,2) - contour(2:end_indx-1, 2)) - (pac(1, 1) - contour(2:end_indx-1, 1))*segment(1,2)) ./ hypot(segment(1, 1), segment(1, 2));
  [junk, indx] = max(dist);
  indx = indx(1);

  if (dist(indx) > thresh)
    indx = indx + 1;
    pac(1:indx, :) = pac_recursive(contour(1:indx, :), thresh);
    pac(indx:end_indx, :) = pac_recursive(contour(indx:end_indx, :), thresh);
  end

  return;
end

function pac = pac_iterative(contour, nsteps, thresh)

  npts = size(contour, 1);
  pac = NaN(size(contour));

  i = 1;
  while (i < npts-nsteps+1)
    pac(i, :) = contour(i, :);
    segment = contour(i + nsteps, :) - contour(i, :);

    dist = abs(segment(1, 1)*(contour(i,2) - contour(i+1:i+nsteps-1, 2)) - (contour(i, 1) - contour(i+1:i+nsteps-1, 1))*segment(1,2)) ./ hypot(segment(1, 1), segment(1, 2));
    [junk, indx] = max(dist);
    indx = indx(1);

    if (dist(indx) > thresh)
      pac(i+indx, :) = contour(i+indx, :);
      i = i + indx;
    else
      i = i + nsteps;
    end
  end

  return;
end
