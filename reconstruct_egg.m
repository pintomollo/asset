function mymovie = reconstruct_egg(mymovie, opts)

  if (nargin == 1 & ischar(mymovie))

    fname = mymovie;
    tmp = dir(fname); 

    for i=1:length(tmp)
      name = tmp(i).name
      load(name);

      mymovie = reconstruct_egg(mymovie, opts);
      save(mymovie.experiment, 'mymovie', 'opts');
    end
    
    return;
  end

  [nframes, imgsize] = size_data(mymovie.dic);
  mymovie = parse_metadata(mymovie, opts);
  real_z = mymovie.metadata.z_position;

  pts = cell(nframes, 1);
  all_coords = zeros(0,3);
  goods = true(nframes, 1);
  sharpness = zeros(nframes, 1);
  all_edges = cell(nframes, 1);

  if (true)
  for i=1:nframes
    tmp_pts = mymovie.dic.eggshell(i).carth;
    if (~isempty(tmp_pts))
      npts = size(tmp_pts, 1);

      img = imnorm(double(load_data(mymovie.dic, i)));
      img = imadm(img, 0, false);
      edges = perpendicular_sampling(img, tmp_pts, opts);
      edges = max(edges(:, 55:75), [], 2);

      sharpness(i) = median(edges);

      tmp_pts = tmp_pts * opts.pixel_size;
      pts{i} = tmp_pts;

      all_edges{i} = edges;
    else
      goods(i) = false;
    end
  end

  save('z_stuff.mat', 'pts', 'all_edges', 'goods');

  %keyboard
  else

  load('z_stuff.mat')
  end
  
  keyboard

  g = true(nframes, 1);
  prev = false(size(g));

  e = [];
  cen = NaN(3, 0);
  a = NaN(3, 0);
  o = NaN(1, 0);

  c = 0;
  while(any(g) & any(xor(prev,g)))
    c = c + 1;
    prev = g;
    [cen(:,c),a(:,c),o(c),g(g),e(c)] = filter_planes(pts(g), all_edges(g), real_z(g));

    if ((c>1 & e(c) > e(c-1)) | (c>2 & (e(c-2) - e(c-1) > e(c - 1) - e(c))))
      c = c-1;
      break;
    end
  end

  mymovie.metadata.center_3d = cen(:,c);
  mymovie.metadata.axes_length_3d = a(:,c);
  mymovie.metadata.orientation_3d = o(1,c);
  
  return;
end

function [fit_center, fit_axes, orient, valids, err] = filter_planes(paths, edges, real_z)

  all_coords = zeros(0,5);
  neggs = length(paths);
  valids = true(neggs,1);

  for i=1:neggs
    if (~isempty(paths{i}))
      all_coords = [all_coords; [paths{i} ones(size(paths{i}, 1), 1)*real_z(i) edges{i} ones(size(edges{i}))*i]];
    else
      valids(i) = false;
    end
  end

  all_coords(:,4) = imnorm(all_coords(:,4));
  val_thresh = graythresh(all_coords(:,4));

  goods = (all_coords(:,4) >= val_thresh);
  goods = filter_goods(goods, 15);

  for i=1:neggs
    current = (all_coords(:,end) == i);
    if ((sum(goods(current)) / sum(current)) < 0.25)
      goods(current) = false;
      valids(i) = false;
    end
  end

  [c,a,o,e] = estimate_axes(all_coords, real_z, goods, valids);
  [m,s] = mymean(e);
  valids = (e < m + 2*s);

  [fit_center,fit_axes,orient,e] = estimate_axes(all_coords, real_z, goods, valids);
  err = mymean(e);

  return;
end

function [fit_center, fit_axes, orient, errs] = estimate_axes(all_coords, real_z, goods, valids)

  neggs = length(real_z);
  currents = ismember(all_coords(:,end), find(valids));
  goods = goods & currents;

  [c, a, orient] = fit_ellipse(all_coords(goods, 1:2));
  new_coords = realign(all_coords(:, 1:2), [0;0], c, orient);
  ell_coords = carth2elliptic(new_coords, [0;0], a, 0);

  ell_coords(:,1) = ell_coords(:,1) - pi;
  [c2,a2,o2] = fit_ellipse(all_coords(:,3),ell_coords(:,2) .* sign(ell_coords(:,1)));

  fit_axes = [a; a2(1)];
  fit_center = [c; c2(1)];

  p0 = [fit_center; fit_axes; orient];
  optims = optimset('Display', 'off');
  lbound = [0;0;-Inf;0;0;0;-2*pi];
  ubound = [Inf(6,1);2*pi];
  bests = lsqcurvefit(@fit_3d_ellipse, p0, all_coords, zeros(neggs, 1), lbound, ubound, optims);

  fit_center = bests(1:3);
  fit_axes = bests(4:6);
  orient = bests(7);

  errs = ellipse_distance(all_coords(goods, :), real_z, neggs, fit_center, fit_axes, orient);

  return;

  function errs = fit_3d_ellipse(p, pts)
    
    errs = ellipse_distance(pts(goods, :), real_z, neggs, p(1:3), p(4:6), p(7));
    errs(~isfinite(errs)) = 0;

    return;
  end
end

function errs = ellipse_distance(coords, z_pos, neggs, center, axes_length, orient)

  if (isempty(coords))
    errs = Inf(neggs,1);

    return;
  end

  valids = unique(coords(:, end)).';
  new_coords = realign(coords(:, 1:2), [0;0], center(1:2, 1), orient);

  z_coef = sqrt(1 - ((z_pos - center(3)).^2 / axes_length(3).^2));
  new_axes = axes_length(1:2, 1) * z_coef;
  errs = NaN(neggs, 1);

  for i=valids
    if (imag(z_coef(i)))
      errs(i) = 1;
    else
      tmp_pts = new_coords(coords(:,end) == i,:);
      ell_pts = carth2elliptic(tmp_pts, [0;0], new_axes(:, i), 0, 'radial');
      errs(i, 1) = mean(abs(ell_pts(:,2) - 1));
    end
  end

  return;
end

function goods = filter_goods(goods, thresh)

  [pos, lengths, vals] = boolean_domains(goods);
  rem = ((lengths < thresh) & ~vals);

  if (any(rem))
    for i=length(rem)-1:-1:2
      if (rem(i))
        lengths(i-1) = lengths(i-1) + lengths(i) + lengths(i+1);
      end
    end
  end

  keep = ~(rem | [false; rem(1:end-1)]);
  pos = pos(keep);
  lengths = lengths(keep);
  vals = vals(keep);

  rem = ((lengths < thresh) & vals);
  vals(rem) = false;
  goods(:) = false;

  for i=1:length(pos)
    if (vals(i))
      goods(pos(i):pos(i)+lengths(i)-1) = true;
    end
  end

  return;
end

function plot_3d_ellipse(pts, z_pos, orient)
  
  if (nargin == 3)
    [centers, axes] = deal(pts, z_pos);

    nslices = 7;
    z_pos = [-nslices:nslices] * (axes(3)) / (nslices+1) ;
    z_coef = sqrt(1 - (z_pos.^2 / axes(3).^2));
    z_pos = z_pos + centers(3);

    for i=1:length(z_pos)
      [x, y] = draw_ellipse(centers(1:2), axes(1:2)*z_coef(i), orient);
      plot3(x,y,z_pos(i)*ones(size(x)), 'r');
    end
  else

    nframes = length(pts);

    for i=1:nframes
      path = pts{i};
      if (isempty(path))
        continue;
      end
      
      [c, a, o] = fit_ellipse(path);
      [x, y] = draw_ellipse(c, a, o);
      plot3(x,y,z_pos(i)*ones(size(x)), 'b');
      plot3(path(:,1),path(:,2),z_pos(i)*ones(size(path(:,1))), 'k');
    end
  end

  return;
end
