function [new_spots] = lidke_fit(img, opts, nimg)

  if (isstruct(img) & isfield(img, 'experiment'))

    mymovie = img;
    nframes = size_data(mymovie.data);

    if (isfield(mymovie.data, 'spots') & ~isempty(mymovie.data.spots))
      spots = mymovie.data.spots;
    else
      spots = get_struct('ruffles', [1,nframes]);
    end

    lastwarn('');

    cpb = ConsoleProgressBar();
    cpb.setLeftMargin(4);   % progress bar left margin
    cpb.setTopMargin(1);    % rows margin
    cpb.setLength(100);      % progress bar length: [.....]
    cpb.setMinimum(0);
    cpb.setMaximum(nframes);

    cpb.setElapsedTimeVisible(1);
    cpb.setRemainedTimeVisible(1);
    cpb.setElapsedTimePosition('left');
    cpb.setRemainedTimePosition('right');

    cpb.start();

    %indxs = randperm(nframes);
    for i=1:nframes
      nimg = i;
      %nimg = indxs(i);
      %nimg = 4
      img = load_data(mymovie.data, nimg);
      pts = lidke_fit(double(img), opts, nimg);

      spots(nimg).cluster = [];
      if (isempty(pts))
        spots(nimg).carth = [];
        spots(nimg).properties = [];
      else
        spots(nimg).carth = pts(:,1:2);
        spots(nimg).properties = pts(:,3:end);
      end

      text = sprintf('Progress: %d/%d', i, nframes);
      cpb.setValue(i);
      cpb.setText(text);
    end

    cpb.stop();

    mymovie.data.spots = spots;
    new_spots = mymovie;
    
    return;
  end

  % Get the maximum spot size (in um), used only to reduce the amount of computation
  spot_max_size = opts.spot_tracking.max_size;

  % If the size of the puxels has been set, we can compute the actual spot size
  if (opts.pixel_size > 0)
    
    % The maximal size of the spot in pixels
    spot_max_size = ceil(spot_max_size / opts.pixel_size);
  else
    
    % Try to compute the pixel size
    opts = set_pixel_size(opts);

    % If it worked, we can now compute the size in pixels
    if (opts.pixel_size > 0)
      spot_max_size = ceil(spot_max_size / opts.pixel_size);
    else
      % Otherwise we take it as such
      spot_max_size = ceil(spot_max_size);
    end
  end
  spot_max_size = ceil(spot_max_size);

  % Get the other parameters used in the detection
  % The first one is used to remove noise from the signal
  noise_thresh = opts.spot_tracking.noise_thresh;

  nmax = opts.spot_tracking.max_particles;
  niter = opts.spot_tracking.max_iterations;
  iter_thresh = opts.spot_tracking.iteration_threshold;
  fit_sep = opts.spot_tracking.fit_intensities;
  fusion_thresh = opts.spot_tracking.fusion_thresh;
  intens_thresh = opts.spot_tracking.intensity_thresh;

%  figure;imagesc(img);

  img2 = imfilter(img, fspecial('gaussian', 7, 0.6), 'symmetric');
%  figure;imagesc(img2);
  atrous = imatrou(img2, spot_max_size);
%  figure;imagesc(atrous);
  params = estimate_noise(img);

  [atAvg, atStd] = localAvgStd2D(atrous, 3*spot_max_size + 1);
  thresh = (atrous > (atAvg + noise_thresh * atStd) & img2 > (params(1) + noise_thresh*params(2)));
%  figure;imagesc(thresh);

  thresh = bwmorph(thresh, 'clean');
  bw = locmax2d(img2 .* thresh, ceil(spot_max_size([1 1])/3));

%  figure;imagesc(bw);
  bw = bwmorph(bw, 'shrink', 1);

%  figure;imagesc(bw);
  [cand_y, cand_x] = find(bw);

  nspots = length(cand_x);
  new_spots = NaN(0, 5);
  indexes = NaN(0, 1);

  % First estimation !
  params(end+1) = sum(sum(img2 - params(1))) / nspots;
  window_size = ceil(spot_max_size/6);
  init_pos = [0, 0, window_size];

  for i=1:nspots
    cs = [cand_x(i), cand_y(i)];
    window = get_window(img, cs, window_size);

    pos = fit_multi(window, init_pos, params, spot_max_size/4, 1, niter, iter_thresh, fit_sep);

    if (~isempty(pos))
      pos(:, 1:2) = bsxfun(@plus, cs, pos(:, 1:2));
      new_spots = [new_spots; pos];
    else
      new_size = ceil(window_size/2);

      window = get_window(img, cs, new_size);
      pos = fit_multi(window, init_pos, params, spot_max_size/4, 1, niter, iter_thresh, fit_sep);

      if (~isempty(pos))
        pos(:, 1:2) = bsxfun(@plus, cs, pos(:, 1:2));
        new_spots = [new_spots; pos];
      end
    end

    indexes = [indexes; ones(size(pos, 1), 1)*i];
  end

%  figure;imagesc(img);
%  hold on;
%  scatter(cand_x, cand_y, 'k')
%  scatter(new_spots(:,1), new_spots(:,2), 'w');

  goods = (new_spots(:,4) > intens_thresh*params(2));
  new_spots = fuse_spots(new_spots(goods, :), [cand_x(indexes(goods)) cand_y(indexes(goods))], fusion_thresh);
  nspots = size(new_spots, 1);

%  scatter(new_spots(:,1), new_spots(:,2), 'b');

  dist = sqrt(bsxfun(@minus, new_spots(:,1), new_spots(:,1).').^2 + bsxfun(@minus, new_spots(:,2), new_spots(:,2).').^2);
  rads = new_spots(:,3);
  rads = bsxfun(@plus, rads, rads.') / 2;
  overlap = (dist < 1.5*rads) & ~eye(nspots);
  alones = ~any(overlap, 2);

  spots = new_spots;
  new_spots = spots(alones, :);
  spots = spots(~alones, :);
  overlap = overlap(~alones, ~alones);
  nspots = size(spots, 1);

  for i=1:nspots
    groups = overlap(:,i);
    groups = any(overlap(:, groups), 2) | groups;

    if (any(groups))
      cands = [spots(groups, 1:3)];
      avg = mymean(cands, 1);

      min_pos = min(bsxfun(@minus, cands(:,1:2), cands(:,3)), [], 1);
      max_pos = max(bsxfun(@plus, cands(:,1:2), cands(:,3)), [], 1);

      window_size = ceil(mean(max_pos - min_pos)/2);
      cs = round(avg(1:2));

      cands(:, 1:2) = bsxfun(@minus, cands(:, 1:2), cs);
      window = get_window(img, cs, window_size);
      pos = fit_multi(window, cands, params, spot_max_size/4, 1, niter, iter_thresh, fit_sep);
  %keyboard
      if (~isempty(pos))
        pos(:, 1:2) = bsxfun(@plus, cs, pos(:, 1:2));
        new_spots = [new_spots; pos];
        overlap(:, groups) = false;
      end
    end
  end
  goods = (new_spots(:,4) > intens_thresh*params(2));
  new_spots = fuse_spots(new_spots(goods, :), [cand_x(indexes(goods)) cand_y(indexes(goods))], fusion_thresh);

%  estimated_img = draw_spots(new_spots, img);
%  scatter(new_spots(:,1), new_spots(:,2), 'm');
%  keyboard

  return;
end

function [fused_spots] = fuse_spots(spots, target, dist_thresh)

  %Fuse spots
  dist = sqrt(bsxfun(@minus, spots(:,1), spots(:,1).').^2 + bsxfun(@minus, spots(:,2), spots(:,2).').^2);
  center_dist = sqrt(bsxfun(@minus, spots(:,1), target(:, 1)).^2 + bsxfun(@minus, spots(:,2), target(:, 2)).^2);

  rads = spots(:,3);
  rads = bsxfun(@plus, rads, rads.') / 2;

  signal = spots(:,3).^2 .* spots(:,4);

  fused = (dist < dist_thresh * rads);
  fused_spots = NaN(0,5);

  for i=1:size(spots, 1)
    groups = fused(:,i);
    groups = any(fused(:, groups), 2);

    if (any(groups))
      if (sum(groups) == 1)
        fused_spots = [fused_spots; spots(groups, :)];
      else

        weights = exp(-center_dist(groups) ./ (2*spots(groups, 3).^2));
        weights = (weights / sum(weights)) .* signal(groups);
        weights = (weights / sum(weights));
        fused_spots = [fused_spots; sum(bsxfun(@times, spots(groups, :), weights), 1)];
      end
      fused(:, groups) = false;
    end
  end

  return;
end

function bests = fit_multi(window, pos, params, spot_max_size, nmax, niter, iter_thresh, fit_sep)

  %%%%% Add a constrain on the IC to stabilize them

  size_window = size(window);
  nwindow = (size_window(1) - 1)/2;
  bests = NaN(0, 4);
  prev_pos = NaN(0,3);
  dist_thresh = nwindow*1.33;
  inner_thresh = nwindow*0.67;

  ninit = size(pos, 1);
  if (ninit > nmax)
    nmax = ninit;
  end

  init_sigma = nwindow/2;
  %iter_thresh = 1e-6;
  init_pos = [0, 0, init_sigma];

  likelihood = NaN(1, nmax+1);

  %% Constraining the solutions
  maxvals = [2 2 1]*spot_max_size;
  minvals = [-2 -2 0.01]*spot_max_size;

  % Computing the likelihood for 0 spots
  uk = params(1);
  if (uk == 0)
    likelihood(1) = 0;
  else
    L0 = window .* (log(uk) - log(window) + 1) - uk;
    L0(window <= 0) = -uk;
    valids = isfinite(L0);
    likelihood(1) = sum(L0(valids));
  end
  best_like = 0;

%  if (ninit > 1)
%  figure;imagesc(window);hold on
%  draw_points(bsxfun(@plus, pos, [1 1 0]*nwindow+1), 'b')
%  end

  for i=ninit:nmax
    for j=1:niter
      [intens, steps, likelihood(i+1)] = step(window, params, pos, nwindow, fit_sep);
      if (all(abs(steps(:)) < iter_thresh))
        break;
      else
        % Undescribed contrains from the original code
        pos = pos - steps/i;

        pos = bsxfun(@min, pos, maxvals);
        pos = bsxfun(@max, pos, minvals);

        if (size(unique(pos(:, 1:2), 'rows'),1) < i)
          break;
        end
      end
%      if (i > 1)
%  draw_points(bsxfun(@plus, pos, [1 1 0]*nwindow+1), 'c')
%  end
    end
%    if (i > 1)
%  draw_points(bsxfun(@plus, pos, [1 1 0]*nwindow+1), 'k')
%  end

    bads = any((pos(:,1:2) < -dist_thresh) | (pos(:,1:2) > dist_thresh), 2);

    p_better = chi2cdf(-2*(likelihood(best_like+1) - likelihood(i+1)), (i - best_like)*4);

    if (p_better > 0.95)
      bests = [pos(~bads, :) intens(~bads)];
      % Rescaling to obtain the amplitude
      bests(:,end) = bests(:,end) ./ (2*pi*(bests(:,3).^2));
      bests = [bests ones(size(bests, 1), 1)*likelihood(i+1)];
      best_like = i;
    else
      break;
    end

    if (all(bads))
      break;
    end

    if (i < nmax)
      if (any(bads))
        if (i==ninit)
          pos = pos(~bads, :);
          nmiss = (i - size(pos, 1));
          inits = repmat(init_pos, nmiss, 1);
          inits(:, 1:2) = inits(:, 1:2) + 2*rand(nmiss, 2)-0.5;
          pos = [pos; inits];
          intens = intens(~bads);
          intens(end+1:end+nmiss, 1) = mean(intens);
        else
          pos = [prev_pos; init_pos];
        end
        prev_pos = pos;
      else
        prev_pos = pos;
      end

      curr_window = window;
      for j=1:i
        curr_window = curr_window - intens(j)*GaussMask2D(pos(j,3), size_window, pos(j,[2 1]), 2);
      end
      curr_window(~isfinite(curr_window)) = -Inf;
      [junk, new] = max(curr_window(:));
      [new_y, new_x] = ind2sub(size_window, new);

      new = [new_x new_y] - nwindow - 1;
      if (any(new < -inner_thresh | new > inner_thresh))
        corrs = (init_sigma/2)*sign(new) + new;
      else
        corrs = (init_sigma/2)*sign(mean(pos(:, 1:2), 1) - new) + new;
      end
      corrs = [corrs init_sigma];

      pos = [pos; corrs];

      if (size(unique(pos(:, 1:2), 'rows'),1) < i+1)
        break;
      end
    end
  end

  return;
end

function window = get_window(img, pos, wsize)
% This function extracts a window of size wsize in img around pos.

  % Get the size of the full image
  [h,w] = size(img);

  % Initialize the window
  window = NaN(2*wsize + 1);

  % Create an index vector to compute which pixels lie outside the image 
  indx = [1:2*wsize+1];

  % First check along the x coordinate
  indxx = indx - wsize - 1 + pos(1);
  % We keep only the valid indexes
  okx = (indxx > 0 & indxx <= w);
  % Then along the y
  indxy = indx - wsize - 1 + pos(2);
  oky = (indxy > 0 & indxy <= h);

  % Assign the pixels from the image to the window. All the outside ones will be NaN
  window(indx(oky), indx(okx)) = img(indxy(oky), indxx(okx));

  return;
end

function [intens, steps, likelihood] = step(img, params, cands, wsize, fit_sep)

  steps = NaN(size(cands));

  % Max value for jumps x, y, o, I, b
  maxjump = [1e0, 1e0, 5e-1, 2e0, 1e2];

  posx = [-wsize:wsize];
  posy = posx.';

  npos = length(posx);
  ncands = size(cands, 1);
  uk = zeros(npos, npos, ncands);

  dukdx = NaN([npos, npos, ncands]);
  dukdy = NaN([npos, npos, ncands]);
  dukdo = NaN([npos, npos, ncands]);
  ddukdx2 = NaN([npos, npos, ncands]);
  ddukdy2 = NaN([npos, npos, ncands]);
  ddukdo2 = NaN([npos, npos, ncands]);

  for i=1:ncands
    % intesity missing, fit it and introduce it after
    norm = (1/(sqrt(2*pi)*cands(i,3)));
    sigma = 1 / (cands(i,3)^2);

    posxp = posx - cands(i,1) + 0.5;
    posxm = posx - cands(i,1) - 0.5;
    posyp = posy - cands(i,2) + 0.5;
    posym = posy - cands(i,2) - 0.5;
    Ex = emissions(cands(i,1), cands(i,3), posx);
    Ey = emissions(cands(i,2), cands(i,3), posy);
    expxp = exp(-(posxp.^2) * sigma * 0.5);
    expxm = exp(-(posxm.^2) * sigma * 0.5);
    expyp = exp(-(posyp.^2) * sigma * 0.5);
    expym = exp(-(posym.^2) * sigma * 0.5);
    ex = (posxm.*expxm - posxp.*expxp);
    ey = (posym.*expym - posyp.*expyp);
    exEy = bsxfun(@times, ex, Ey);
    eyEx = bsxfun(@times, ey, Ex);

    dukdx(:,:,i) = norm*bsxfun(@times, Ey, expxm - expxp);
    dukdy(:,:,i) = norm*bsxfun(@times, expym - expyp, Ex);

    ddukdx2(:,:,i) = norm*sigma*exEy;
    ddukdy2(:,:,i) = norm*sigma*eyEx;

    dukdo(:,:,i) = (ddukdx2(:,:,i) + ddukdy2(:,:,i)) * cands(i,3);
    ddukdo2(:,:,i) = norm*(sigma^2)*(bsxfun(@times, ((posxm.^3).*expxm - (posxp.^3).*expxp), Ey) + bsxfun(@times, ((posym.^3).*expym - (posyp.^3).*expyp), Ex)) - 2*dukdo(:,:,i)/cands(i,3);

    uk(:,:,i) = bsxfun(@times, Ey, Ex);
  end
  
  intens = NaN(ncands, 1);

  switch (fit_sep)
    case 'separate'
      % Full fitting
      alls = reshape(uk, npos^2, []);
      valids = isfinite(img);
      intens = (alls(valids, :) \ (img(valids) - params(2)));
      intens_thresh = max(intens) / 5;
      intens(intens < intens_thresh) = intens_thresh;
      uk = reshape(sum(bsxfun(@times, alls, intens.'), 2), npos, npos) + params(2);
    case 'together'
      % Single fit
      uk = sum(uk, 3);
      valids = isfinite(img);
      intens(:) = abs(uk(valids) \ (img(valids) - params(1)));
      uk = sum(uk, 3)*intens(1) + params(1);
    otherwise
      intens(:) = params(end);
      uk = sum(uk, 3)*intens(1) + params(1);
  end

  % Undescribed contrains from the original code
  cf = ((img ./ uk) - 1);
  df = (img ./ (uk.^2));

  cf(uk < 10e-3) = 0;
  df(uk < 10e-3) = 0;

  cf = min(cf, 10e4);
  df = min(df, 10e4);
  % Undescribed contrains from the original code

  for i=1:ncands
    xnum = dukdx(:,:,i) .* cf;
    xdenom = (ddukdx2(:,:,i) .* cf - ((dukdx(:,:,i).^2) .* df * intens(i)));
    stepx = sum(xnum(isfinite(xnum))) / sum(xdenom(isfinite(xdenom)));

    ynum = dukdy(:,:,i) .* cf;
    ydenom = (ddukdy2(:,:,i) .* cf - ((dukdy(:,:,i).^2) .* df * intens(i)));
    stepy = sum(ynum(isfinite(ynum))) / sum(ydenom(isfinite(ydenom)));

    onum = dukdo(:,:,i) .* cf;
    odenom = (ddukdo2(:,:,i) .* cf - ((dukdo(:,:,i).^2) .* df * intens(i)));
    stepo = sum(onum(isfinite(onum))) / sum(odenom(isfinite(odenom)));

    steps(i,:) = 0.5*[stepx, stepy, stepo];
  end

  % Bounding the moves (undescribed)
  steps = bsxfun(@min, steps, maxjump(1:3));

  likelihood = img .* (log(uk) - log(img) + 1) - uk;
  % Undescribed contrains from the original code
  likelihood(uk <= 0) = 0;
  likelihood(uk > 0 & img <= 0) = -uk(uk > 0 & img <= 0);
  % Undescribed contrains from the original code
  valids = isfinite(likelihood);

  steps(isnan(steps)) = 0;
  steps = steps;

  likelihood = sum(likelihood(valids));

  return
end

function intens = emissions(p0, sigma, p)

  norm = 1/(sqrt(2)*sigma);
  intens = 0.5 * (erf((p - p0 + 0.5)*norm) - erf((p - p0 - 0.5)*norm));

  return;
end

function draw_points(pts, color, min_intens)
  
  if (nargin == 1)
    color = 'b';
%    min_intens = min(pts(:,4));
  elseif (nargin == 2)
%    min_intens = min(pts(:,4));
  end

%  intens = pts(:,4) / (2*min_intens);

  for i=1:size(pts, 1)
    %rectangle('Position', [pts(i,1)-pts(i, 3), pts(i,2) - pts(i,3), max(2*pts(i,[3 3]), [1e-1 1e-1])], 'Curvature', [1 1], 'EdgeColor', color, 'LineWidth', intens(i));
    rectangle('Position', [pts(i,1)-pts(i, 3), pts(i,2) - pts(i,3), max(2*pts(i,[3 3]), [1e-1 1e-1])], 'Curvature', [1 1], 'EdgeColor', color);
  end

  return;
end
