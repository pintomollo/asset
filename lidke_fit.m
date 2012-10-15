function lidke_fit(img)

  nspots = 200;
  noise_thresh = 0;
  sigma = 4;
  spot_max_size = 3*sigma;
  nmax = 5;

  niter = 20;
  iter_thresh = 1e-10;

  params = [500 10 sigma];


  if (true)
    if (nargin == 0)
      img = zeros(600, 400);
      size_img = size(img);
      myspots = bsxfun(@times, rand(nspots, 2), size_img([2 1]));
      myspots = [myspots (((rand(nspots, 1)-0.5)/2)+1)*params(3)];

      for i=1:nspots
        img = img  + GaussMask2D(myspots(i, 3), size_img, myspots(i, [2 1]), 2, 1) * params(1);
      end
      img = img + params(2);
      %img = img + poissrnd(1, size(img));
      %img = img + randn(size(img))*3;
    end

    img2 = imfilter(img, fspecial('gaussian', 7, 0.6), 'symmetric');
    atrous = imatrou(img2, noise_thresh, spot_max_size);
    thresh = mean(atrous(:)) + noise_thresh * std(atrous(:));
    bw = imfill(atrous > thresh, 'holes');
    bw = bwmorph(bw, 'shrink', Inf);
    [cand_y, cand_x] = find(bw);
  else
    img = zeros(600, 400);
    size_img = size(img);
    myspots = [60 30; 90 20];
    for i=1:size(myspots, 1)
      img = img  + GaussMask2D(sigma*2.5, size_img, myspots(i, [2 1]), 2, 1);
    end
    img = img*params(1) + params(2);
    %img = img + poissrnd(1, size(img));

    [Y,X] = meshgrid([-sigma:2:sigma] + myspots(1,2),[-sigma:3:sigma] + myspots(1,1));
    cand_y = Y(:);
    cand_x = X(:);
    %cand_x = myspots(:,1);
    %cand_y = myspots(:,2);
  end

  nspots = length(cand_x);

  new_spots = NaN(0, 3);
  
  %figure;

  for i=1:nspots
    cs = [cand_x(i), cand_y(i)];
    window = get_window(img, cs, spot_max_size);

    pos = fit_multi(window, params, spot_max_size, nmax, niter, false);

    if (~isempty(pos))
      new_spots = [new_spots; bsxfun(@plus, [cs 0], pos)];
    else
      window = get_window(img, cs, spot_max_size/2);
      pos = fit_multi(window, params, spot_max_size/2, nmax, niter, false);

      if (~isempty(pos))
        new_spots = [new_spots; bsxfun(@plus, [cs 0], pos)];
      end
    end
  end

  %Fuse spots
  dist = sqrt(bsxfun(@minus, new_spots(:,1), new_spots(:,1).').^2 + bsxfun(@minus, new_spots(:,2), new_spots(:,2).').^2);
  rads = params(3)/2;

  fused = dist < rads;
  fused_spots = NaN(0,3);
  for i=1:size(new_spots, 1)
    groups = fused(:,i);

    if (any(groups))
      fused_spots = [fused_spots; mymean(new_spots(groups, :), 1)];
      fused(:, groups) = false;
    end
  end

  figure;imagesc(img);hold on;
  scatter(myspots(:,1), myspots(:,2), 'w')
  scatter(cand_x, cand_y, 'y')
  scatter(new_spots(:,1), new_spots(:,2), 'r');
  scatter(fused_spots(:,1), fused_spots(:,2), 'k');

  return;
end

function bests = fit_multi(window, params, spot_max_size, nmax, niter, show_it)

  p_thresh = 0.5;
  iter_thresh = 1e-6;
  pos = [0, 0, params(3)];
  size_window = size(window);
  max_like = -Inf;
  bests = [];
  prev_pos = NaN(0,3);
  dist_thresh = spot_max_size*1.33;
  inner_thresh = spot_max_size*0.67;

  %% Constraining the solutions
  maxvals = [2 2 0.75]*spot_max_size;
  minvals = [-2 -2 0.01]*spot_max_size;

  for i=1:nmax
    if (show_it)
      hold off;imagesc(window);hold on;
      %scatter(pos(:,1)+spot_max_size+1, pos(:,2)+spot_max_size+1, 'y');
      draw_points([pos(:,1)+spot_max_size+1, pos(:,2)+spot_max_size+1 pos(:,3)], 'y');
    end

    for j=1:niter
      [params, steps, likelihood] = step(window, params, pos, spot_max_size);
      if (all(abs(steps(:)) < iter_thresh))
        break;
      else
        %pos = pos - steps;
        % Undescribed contrains from the original code
        pos = pos - steps/i;
        
        pos = bsxfun(@min, pos, maxvals);
        pos = bsxfun(@max, pos, minvals);

        if (show_it)
          %scatter(pos(:,1)+spot_max_size+1, pos(:,2)+spot_max_size+1, 'c');
          draw_points([pos(:,1)+spot_max_size+1, pos(:,2)+spot_max_size+1 pos(:,3)], 'c');
          %pause(0.5)
        end
      end
    end

    bads = any((pos(:,1:2) < -dist_thresh) | (pos(:,1:2) > dist_thresh), 2);

    if (likelihood > max_like & likelihood > p_thresh)
      bests = pos(~bads, :);
      max_like = likelihood;

      if (show_it)
        draw_points([pos(:,1)+spot_max_size+1, pos(:,2)+spot_max_size+1 pos(:,3)], 'k');
        %scatter(pos(:,1)+spot_max_size+1, pos(:,2)+spot_max_size+1, 'k');
        title([num2str(likelihood) ' (' num2str(i) ') !!'])
        pause(1)
      end

      break;
    end
    if (show_it)
      %scatter(pos(:,1)+spot_max_size+1, pos(:,2)+spot_max_size+1, 'k');
      draw_points([pos(:,1)+spot_max_size+1, pos(:,2)+spot_max_size+1 pos(:,3)], 'k');
      title([num2str(likelihood) ' (' num2str(i) ')'])
    end

    if (all(bads))
      pos = [prev_pos; [0,0, params(3)]];
      if (i>1 & ~any(pos(:)))
        pos(1:end-1, 1:2) = (rand(i-1, 2) - 0.5)*spot_max_size;
      end
      prev_pos = pos;
    else
      prev_pos = pos;
    end

    curr_window = window;
    for j=1:i
      curr_window = curr_window - params(1)*GaussMask2D(pos(j,3), size_window, pos(j,[2 1]), 2);
    end
    curr_window(~isfinite(curr_window)) = -Inf;
    [junk, new] = max(curr_window(:));
    [new_y, new_x] = ind2sub(size_window, new);

    new = [new_x new_y] - spot_max_size - 1;
    if (any(new < -inner_thresh | new > inner_thresh))
      corrs = (params(3)/2)*sign(new) + new;
    else
      corrs = (params(3)/2)*sign(mean(pos(:, 1:2), 1) - new) + new;
    end
    corrs = [corrs params(3)];

    if (show_it)
      scatter(new(1)+spot_max_size+1, new(2)+spot_max_size+1, 'm');
      scatter(corrs(1)+spot_max_size+1, corrs(2)+spot_max_size+1, 'r');
      pause(1)
    end

    pos = [pos; corrs];
    %keyboard
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

function [params, steps, likelihood] = step(img, params, cands, wsize)

  steps = NaN(size(cands));

  % Max value for jumps x, y, o, I, b
  maxjump = [1e0, 1e0, 5e-1, 2e0, 1e2];

  posx = [-wsize:wsize];
  posy = posx.';

  npos = length(posx);
  ncands = size(cands, 1);
  uk = zeros(npos, npos);

  dukdx = NaN([npos, npos, ncands]);
  dukdy = NaN([npos, npos, ncands]);
  dukdo = NaN([npos, npos, ncands]);
  dukdI = NaN([npos, npos, ncands]);
  dukdb = ones(npos, npos);
  ddukdx2 = NaN([npos, npos, ncands]);
  ddukdy2 = NaN([npos, npos, ncands]);
  ddukdI2 = zeros(npos, npos);
  ddukdb2 = zeros(npos, npos);
  ddukdo2 = NaN([npos, npos, ncands]);

  for i=1:ncands
    norm = (params(1)/(sqrt(2*pi)*cands(i,3)));
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
    dukdI(:,:,i) = bsxfun(@times, Ey, Ex);

    ddukdx2(:,:,i) = norm*sigma*exEy;
    ddukdy2(:,:,i) = norm*sigma*eyEx;

    dukdo(:,:,i) = (ddukdx2(:,:,i) + ddukdy2(:,:,i)) * cands(i,3);
    ddukdo2(:,:,i) = norm*(sigma^2)*(bsxfun(@times, ((posxm.^3).*expxm - (posxp.^3).*expxp), Ey) + bsxfun(@times, ((posym.^3).*expym - (posyp.^3).*expyp), Ex)) - 2*dukdo(:,:,i)/cands(i,3);
    %2*(norm/(sqrt(2*pi)))*cands(i,3)*sigma*bsxfun(@times, ex, ey);

    %norm*(sigma^2)*(sigma*(bsxfun(@times, ((posxm.^3).*expxm - (posxp.^3).*expxp), Ey) + bsxfun(@times, ((posym.^3).*expym - (posyp.^3).*expyp), Ex)) - ...
    %                4*(exEy + eyEx) + 2*norm*bsxfun(@times, ey, ex));
  end

  uk = sum(dukdI, 3)*params(1) + params(2);

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
    xdenom = (ddukdx2(:,:,i) .* cf - ((dukdx(:,:,i).^2) .* df));
    stepx = sum(xnum(isfinite(xnum))) / sum(xdenom(isfinite(xdenom)));

    ynum = dukdy(:,:,i) .* cf;
    ydenom = (ddukdy2(:,:,i) .* cf - ((dukdy(:,:,i).^2) .* df));
    stepy = sum(ynum(isfinite(ynum))) / sum(ydenom(isfinite(ydenom)));

    inum = dukdI(:,:,i) .* cf;
    idenom = - ((dukdI(:,:,i).^2) .* df);
    stepi = sum(inum(isfinite(inum))) / sum(idenom(isfinite(idenom)));

    onum = dukdo(:,:,i) .* cf;
    odenom = (ddukdo2(:,:,i) .* cf - ((dukdo(:,:,i).^2) .* df));
    stepo = sum(onum(isfinite(onum))) / sum(odenom(isfinite(odenom)));

    %steps(i,:) = [stepx, stepy, stepo, stepi];
    %steps(i,:) = 0.5*[stepx, stepy, 0];
    steps(i,:) = 0.5*[stepx, stepy, stepo];
  end

  % Bounding the moves (undescribed)
  steps = bsxfun(@min, steps, maxjump(1:3));

  bnum = cf;
  bdenom =  -df;
  stepb = sum(inum(isfinite(inum))) / sum(idenom(isfinite(idenom)));
  stepb = min(stepb, maxjump(end));

  likelihood = img .* (log(uk) - log(img) + 1) - uk;
  % Undescribed contrains from the original code
  likelihood(uk <= 0) = 0;
  likelihood(uk > 0 & img <= 0) = -uk(uk > 0 & img <= 0);
  % Undescribed contrains from the original code
  valids = isfinite(likelihood);

  %noise = 0;
  %best_like = noise*sum(valids(:)) - sum((uk(valids) + noise).*log(1+noise./uk(valids)));
  best_like = 0;

  steps(isnan(steps)) = 0;
  steps = steps;

  likelihood = exp(2*(sum(likelihood(valids)) - best_like));

  return
end

function intens = emissions(p0, sigma, p)

  norm = 1/(sqrt(2)*sigma);
  intens = 0.5 * (erf((p - p0 + 0.5)*norm) - erf((p - p0 - 0.5)*norm));

  return;
end

function draw_points(pts, color)
  
  if (nargin == 1)
    color = 'b';
  end

  for i=1:size(pts, 1)
    rectangle('Position', [pts(i,1)-pts(i, 3), pts(i,2) - pts(i,3), max(2*pts(i,[3 3]), [1e-1 1e-1])], 'Curvature', [1 1], 'EdgeColor', color);
  end

  return;
end
