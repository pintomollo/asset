function [pac, hits] = impac(contour, nsteps, thresh, method)

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

  if (all(contour(1,:) ~= contour(end, :)))
    contour = contour([1:end 1], :);
  end
  npts = size(contour, 1);

  if (isempty(thresh))
    % Rough incircle radius estimation (from the bounding box)
    thresh = 0.1 * min((max(contour, [], 1) - min(contour, [], 1))) / 2;

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
      indx = ceil(size(contour) / 2);
      pac(1:indx, :) = pac_recursive(contour(1:indx, :), thresh);
      pac(indx:end, :) = pac_recursive(contour(indx:end, :), thresh);
    case 'iterative'
      pac = pac_iterative(contour, nsteps, thresh);
  end

  pac = pac([1:end-1], :);
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

  if (end_indx == 1)
    return;
  end

  pac([1 end_indx], :) = contour([1 end_indx], :);

  dist = abs(segment(1, 1)*(pac(1,2) - contour(2:end_indx-1, 2)) - (pac(1, 1) - contour(2:end_indx-1, 1))*segment(1,2)) ./ hypot(segment(1, 1), segment(1, 2));
  [~, indx] = max(dist);
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
    [~, indx] = max(dist);
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
