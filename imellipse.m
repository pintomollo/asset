function ellipse = imellipse(img, wsize, thresh, angle_thresh, radius_thresh, shortcut)

  size_img = size(img);

  if (nargin == 1)
    thresh = 0;
    wsize = 2;
    angle_thresh = pi/8;
  elseif (nargin == 2)
    thresh = 0;
  end

  %keyboard

  alpha = pi/4;
  precision_thresh = 1e-2;
  voting_thresh = 50;

  [edges, direct] = imadm(img, thresh, true);

  binary_edges = double(edges > 0);

  neigh_filter = [0 1 0; 1 0 1; 0 1 0];
  neighborhood = imfilter(binary_edges, neigh_filter, 'symmetric');
  binary_edges(neighborhood == 0 | neighborhood > 2) = 0;

  filter_size = 2*wsize + 1;
  positives = ones(filter_size);
  negatives = -positives;
  nulls = zeros(filter_size);

  hfilt = nulls;
  hfilt(1:wsize, :) = negatives(1:wsize, :);
  hfilt(wsize+2:end, :) = positives(wsize+2:end, :);
  vfilt = -hfilt.';

  ndfilt = triu(positives, 1) + tril(negatives, -1);
  pdfilt = fliplr(ndfilt);

  himg = imfilter(edges, hfilt, 'symmetric');
  vimg = imfilter(edges, vfilt, 'symmetric');
  pdimg = imfilter(edges, pdfilt, 'symmetric');
  ndimg = imfilter(edges, ndfilt, 'symmetric');

  convexity = cat(3, ndimg, vimg, pdimg, himg);
  angles = round(direct / (pi/4));
  angles(angles > 4) = 4 - angles(angles > 4);
  angles(angles == 0) = 4;

  inds = sub2ind(size(convexity), ...
                 repmat([1:size(img, 1)].', 1, size(img, 2)), ...
                 repmat(1:size(img, 2), size(img, 1), 1), ...
                 abs(angles));
  convexity = convexity(inds);
  % Add a smoothing term to avoid having both sides of the line with different convexities
  convexity = imfilter(convexity, positives, 'symmetric');

  % We could even filter out the small convexities (might lead to problems ?)

  convex_angle = direct + (pi/2) * sign(angles) .* sign(convexity);
  convex_angle(convex_angle < 0) = convex_angle(convex_angle < 0) + 2*pi;
  convex_angle(convex_angle > 2*pi) = convex_angle(convex_angle > 2*pi) - 2*pi;
  convex_angle(convexity == 0 | binary_edges == 0) = NaN;

  clearvars himg vimg pdimg ndimg convexity angles inds neighborhood;
  
  %figure;imagesc(convex_angle);

  % Compute the actual angle !
  img_size = size(img); 

  indexes = NaN(0, 7);
  cands = cell(numel(img), 2);
  %count = 0;

  done = false(size(img));
  %found = done;
  first = true;
  %keyboard

  for i = 1:img_size(1)
    if (~any(binary_edges(i, :)))
      continue;
    end

    for j = 1:img_size(2)
      if (done(i,j) & shortcut)
      %  count = count + 1;
        continue;
      end

      % Here we could filter out previously found pts as in the article
      if (binary_edges(i, j) == 1 & ~isnan(convex_angle(i, j)))
          %if (i > 280 & first)
          %  first = false;
          %  keyboard;
          %end

          linear_scan(i, j, -alpha);
          linear_scan(i, j, alpha);
      end
    end
  end

  %keyboard

  window = nulls;
  window(:, wsize+1:end) = 1;
  [tmp_x, tmp_y] = find(window);
  window = (tmp_x - (wsize + 1)) + (tmp_y - (wsize + 1))*img_size(1);
  %window = find(window) - (wsize^2 + (wsize+1)^2);

  indexes = find(done);
  %cands = cands(indexes, :);
  angles = convex_angle(indexes);

  npts = length(indexes);
  pairs = zeros(0, 4);

  for i = 1:npts-1
    bounds = window + indexes(i);
    sub_list = indexes(i+1:end);
    neighbors = ismember(sub_list, bounds);
    sub_list = sub_list(neighbors, :);
    sub_angles = angles(i+1:end);
    sub_angles = sub_angles(neighbors);
    curr_angle = angles(i) + [-alpha alpha];

    angle_dist = abs(bsxfun(@minus, sub_angles(:) - alpha, curr_angle));
    angle_dist = min(angle_dist, abs(angle_dist - 2*pi));
    perpend = (abs(angle_dist - pi/2) < precision_thresh);

    [p, q] = find(perpend);
    cache = ones(size(p));
    pairs = [pairs; sub_list(p) cache indexes(i)*cache q];

    angle_dist = abs(bsxfun(@minus, sub_angles(:) + alpha, curr_angle));
    angle_dist = min(angle_dist, abs(angle_dist - 2*pi));
    perpend = (abs(angle_dist - pi/2) < precision_thresh);

    [p, q] = find(perpend);
    cache = ones(size(p));
    pairs = [pairs; sub_list(p) 2*cache indexes(i)*cache q];

    %x_bound = indexes(i, 1) + [0 wsize];
    %x_bound = indexes(i, 1) + [0 wsize];
    %y_bound = indexes(i, 2) + [-wsize wsize];

    %sub_list = indexes(i+1:end);

    %neighbors = (sub_list(:, 1) >= x_bound(1) & sub_list(:, 1) <= x_bound(2)) & (sub_list(:, 2) >= y_bound(1) & sub_list(:, 2) <= y_bound(2));

    %angle_dist = abs(bsxfun(@minus, sub_list(:, 4), indexes(i, 4:5)));
    %angle_dist = min(angle_dist, abs(angle_dist - 2*pi));
    %perpend = (abs(angle_dist - pi/2) < precision_thresh);

    %[p, q] = find(perpend);
    %cache = ones(size(p));
    %pairs = [pairs; sub_list(p, end) cache i*cache q];

    %angle_dist = abs(bsxfun(@minus, sub_list(:, 5), indexes(i, 4:5)));
    %angle_dist = min(angle_dist, abs(angle_dist - 2*pi));
    %perpend = (abs(angle_dist - pi/2) < precision_thresh);

    %[p, q] = find(perpend);
    %cache = ones(size(p));
    %pairs = [pairs; [sub_list(p, end) 2*cache i*cache q]];
  end

  npairs = size(pairs, 1);
  votes = NaN(npairs, 5);

  accum = zeros(img_size);

  for i = 1:npairs
    %indx1 = pairs(i, 1);
    %[tmp_x, tmp_y] = ind2sub(img_size, indx1);
    %pt1 = [tmp_x, tmp_y, convex_angle(indx1)];
    %cands1 = cands{indx1, pairs(i, 2)};

    %indx2 = pairs(i, 3);
    %[tmp_x, tmp_y] = ind2sub(img_size, indx2);
    %pt2 = [tmp_x, tmp_y, convex_angle(indx2)];
    %cands2 = cands{indx2, pairs(i, 4)};

    pt1 = get_pts(pairs(i, 1));
    pt2 = get_pts(pairs(i, 3));
    cands1 = get_pts(cands{pairs(i, 1), pairs(i, 2)});
    cands2 = get_pts(cands{pairs(i, 3), pairs(i, 4)});

    sigma_adapted = wsize + 0.5;
    %sigma_adapted = wsize + 1/mean(indexes([indx1 indx2], 6));

    [tmp_votes] = compute_estimation(pt1, cands1, pt2, cands2);
    if (~isempty(tmp_votes))
      votes(i, :) = tmp_votes;
      %'Need to move out the Gaussian, the paste it at the correct position'
      %for j = 1:size(tmp_votes, 1)
      %  accum = accum + GaussMask2D(sigma_adapted, img_size, tmp_votes(j, [2 1]), 2, 1);
      %end
    end
  end

  %votes(:, end) = align_orientations(votes(:, end));
  %figure;imagesc(accum)

  %keyboard

  return;

  function new_pts = get_pts(pt_indx)
    
    pt_indx = pt_indx(:);
    [tmp1, tmp2] = ind2sub(img_size, pt_indx);
    new_pts = [tmp1, tmp2, convex_angle(pt_indx)];

    return;
  end

  function linear_scan(line_x, line_y, line_angle)
    pt_angle = convex_angle(line_x, line_y);

    % Need to invert the y-axis for the basic geometry to work in the usual coordinate system
    [line_diag] = get_indexes(line_y, line_x, -(pt_angle + line_angle), img_size, radius_thresh);
    
    diag_angles = convex_angle(line_diag);
    diag_good = get_candidates(diag_angles, pi + pt_angle + 2*line_angle, angle_thresh);

    if (false)
      tmp_img = img;
      tmp_img(line_diag) = Inf;
      tmp_img(line_diag(diag_good)) = -Inf;
      imagesc(tmp_img);
      drawnow;
    end

    if (any(diag_good))
      done(line_x, line_y) = true;

      %diag_angles = diag_angles(diag_good);
      line_diag = line_diag(diag_good);
      
      %pos_indx = sub2ind(img_size, line_x, line_y);
      pos_indx = line_x + (line_y-1)*img_size(1);

      if (line_angle < 0)
        cands{pos_indx, 1} = [cands{pos_indx, 1} line_diag];
      else
        cands{pos_indx, 2} = [cands{pos_indx, 2} line_diag];
      end

      already_done = done(line_diag);
      done(line_diag) = true;
      %found(line_diag) = true;

      for tmp_i = 1:length(already_done)

        if (line_angle > 0)
          cands{line_diag(tmp_i), 1} = [cands{line_diag(tmp_i), 1} pos_indx];
        else
          cands{line_diag(tmp_i), 2} = [cands{line_diag(tmp_i), 2} pos_indx];
        end
        
        if (~already_done(tmp_i))
          [pts_x, pts_y] = ind2sub(img_size, line_diag(tmp_i));
          linear_scan(pts_x, pts_y, line_angle);
        end
      end
    end

    return;
  end
end

function [votes] = compute_estimation(pt1, cands1, pt2, cands2)

  ncands1 = size(cands1, 1);
  ncands2 = size(cands2, 1);

  if (ncands1 == 0 | ncands2 == 0)
    votes = [];
    return;
  end

  cands1 = [cands1(:, 2) cands1(:,1) (pi/2)-cands1(:,3)];
  cands2 = [cands2(:, 2) cands2(:,1) (pi/2)-cands2(:,3)];


  pt1 = [pt1(2) pt1(1) (pi/2 - pt1(3))];
  pt2 = [pt2(2) pt2(1) (pi/2 - pt2(3))];

  if (true)

  i = randi(ncands1);
  j = randi(ncands2);

  params1 = triangulate_points(pt1, cands1(i,:));
  params2 = triangulate_points(pt2, cands2(j,:));

  votes = infer_params(params1, params2, [pt1(1:2); pt2(1:2); cands1(i, 1:2); cands2(j, 1:2)]);

  else

  votes = NaN(ncands1*ncands2, 5);
  params1 = zeros(ncands1, 4);
  params2 = zeros(ncands2, 4);

  for i = 1:ncands1
    params1(i, :) = triangulate_points(pt1, cands1(i, :));
  end
  for i = 1:ncands2
    params2(i, :) = triangulate_points(pt2, cands2(i, :));
  end

  for i = 1:ncands1
    for j = 1:ncands2
      votes((i-1)*ncands2 + j, :) = infer_params(params1(i, :), params2(j, :), [pt1(1:2); pt2(1:2); cands1(i, 1:2); cands2(j, 1:2)]);
    end
  end

  end

  return;
end

function params = infer_params(triangle1, triangle2, pts)

  x1 = triangle1(1);
  y1 = triangle1(2);
  q1 = triangle1(3);
  q2 = triangle1(4);

  x2 = triangle2(1);
  y2 = triangle2(2);
  q3 = triangle2(3);
  q4 = triangle2(4);

  a0 = (y2 - q4*x2 - y1 + q2*x1) / (q2 - q4);
  b0 = q2 * (a0 - x1) + y1;
  K = abs(sqrt(1 + ((q1*q2 + 1)*(q3 + q4) - (q3*q4 + 1)*(q1 + q2)) / (q1*q2 - q3*q4)));
  N = abs(sqrt(-((q1 - K)*(q2 - K))/((1 + q1*K)*(1 + q2*K))));

  x0 = ((pts(:,1) - a0) / sqrt(K^2 + 1)) + ((pts(:,2) - b0)*K / sqrt(K^2 + 1));
  y0 = ((pts(:,1) - a0)*K / sqrt(K^2 + 1)) + ((pts(:,2) - b0) / sqrt(K^2 + 1));
  ax = sqrt(((x0.^2)*(N^2) + y0.^2) / (N^2*(1+K^2)));
  ax = mean(abs(ax(:)));

  rho = atan(K);
  a = ax / cos(rho);
  b = N * a;

  params = [a0 b0 a b rho];

  return;
end

function params = triangulate_points(pt1, pt2)

  x1 = pt1(1); 
  y1 = pt1(2); 
  a1 = tan(pt1(3)); 
  x2 = pt2(1); 
  y2 = pt2(2); 
  a2 = tan(pt2(3)); 

  xt = (y2 - a2 * x2 - y1 + a1 * x1) / (a1 - a2);
  yt = a1 * (xt - x1) + y1;

  xm = (x1 + x2) / 2;
  ym = (y1 + y2) / 2;

  q1 = (y2 - y1) / (x2 - x1);
  q2 = (yt - ym) / (xt - xm);

  params = [xm, ym, q1, q2];

  return;
end

function [goods] = get_candidates(angles, objective, thresh)

  % Need to avoid having to check angles
  angles(angles < 0) = angles(angles < 0) + 2*pi;
  angles(angles > 2*pi) = angles(angles > 2*pi) - 2*pi;

  if (objective < 0)
    objective = objective + 2*pi;
  elseif (objective > 2*pi)
    objective = objective - 2*pi;
  end

  % Maybe remove the modulo, might be faster
  %angles = mod(angles, 2*pi); 
  %objective = mod(objective, 2*pi);

  angle_diff = abs(angles - objective);
  goods = (angle_diff < thresh | abs(angle_diff - 2*pi) < thresh);

  return;
end

function [diag] = get_indexes(x, y, angle, maxs, thresh)

  range = [cos(angle) * thresh; sin(angle) * thresh]; 
  range = [floor(range(:, 1)) ceil(range(:, 2))];

  dl = abs(diff(range, 1, 2));
  if (dl(1) < dl(2))
    if (range(2,1) < range(2,2))
      ptsy = [range(2, 1) : range(2, 2)];
    else
      ptsy = [range(2, 2) : range(2, 1)];
    end

    ptsx = ptsy / tan(angle);
    ptsx = [floor(ptsx) ceil(ptsx)];

    ptsy = [ptsy ptsy];
  else
    if (range(1,1) < range(1,2))
      ptsx = [range(1, 1) : range(1, 2)];
    else
      ptsx = [range(1, 2) : range(1, 1)];
    end
  
    ptsy = ptsx * tan(angle);
    ptsy = [floor(ptsy) ceil(ptsy)];

    ptsx = [ptsx ptsx];
  end

  ptsx = ptsx + x;
  ptsy = ptsy + y;

  %'try to avoid this test'
  valids = (ptsx > 0 & ptsx < maxs(2) & ptsy > 0 & ptsy < maxs(1));
  ptsx = ptsx(valids);
  ptsy = ptsy(valids);

  diag = ptsy + (ptsx-1)*maxs(1);
  %diag = sub2ind(maxs, ptsy, ptsx);

  return;
end
