function [mymovie, updated] = reconstruct_egg(mymovie, opts)

  updated = true;
  if (nargin == 1 & ischar(mymovie))

    fname = mymovie;
    tmp = dir(fname); 
    updated = false;

    for i=1:length(tmp)
      name = tmp(i).name
      load(name);

      [mymovie, curr_up] = reconstruct_egg(mymovie, opts);

      if (cur_up)
        save(mymovie.experiment, 'mymovie', 'opts');
        updated = true;
      end
    end
    
    return;
  end

  if ((isfield(mymovie, 'metadata') & isfield(mymovie.metadata, 'relative_z') & ~isempty(mymovie.metadata.relative_z)) & ~opts.recompute)
    updated = false;

    return;
  end

  [nframes, imgsize] = size_data(mymovie.dic);
  mymovie = parse_metadata(mymovie, opts);
  real_z = mymovie.metadata.z_position;

  pts = cell(nframes, 1);
  all_coords = zeros(0,3);
  sharpness = zeros(nframes, 1);
  all_edges = cell(nframes, 1);

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
    end
  end

  if (isempty(mymovie.metadata.plane_index))
    mymovie.metadata.plane_index = ones(1, nframes);
  end
  if (isempty(mymovie.metadata.frame_index))
    mymovie.metadata.frame_index = [1:nframes];
  end

  valids = true(1, nframes);

  for i=1:length(pts)
    if (~isempty(pts{i}))
      if (~isempty(real_z))
        pts{i} = [pts{i} ones(size(pts{i}, 1), 1)*real_z(1, i) ones(size(all_edges{i}))*i all_edges{i}];
      else
        pts{i} = [pts{i} NaN(size(pts{i}, 1), 1) ones(size(all_edges{i}))*i all_edges{i}];
      end
    else
      valids(i) = false;
    end
  end

  frames = mymovie.metadata.frame_index(1,:);
  [time_points, junk, indexes] = unique(frames(~isnan(frames)));
  ntimes = length(time_points);
  centers = NaN(3, ntimes);
  axes_length = NaN(3, ntimes);
  orient = NaN(1, ntimes);
  all_pts = cell(nframes, 1);

  for i=1:ntimes
    currents = (frames == time_points(i) & valids);
    centers(1:2, i) = mean(mymovie.dic.centers(:, currents) * opts.pixel_size, 2);
    axes_length(1:2, i) = mean(mymovie.dic.axes_length(:, currents) * opts.pixel_size, 2);
    orient(1, i) = mean(align_orientations(mymovie.dic.orientations(1, currents)), 2);

    if (sum(currents) > 3)
      [centers(3, i), axes_length(:, i), all_pts(currents)] = fit_absolute_ellipse(pts(currents), centers(1:2, i), axes_length(1:2, i), orient(1, i));

    elseif (sum(currents) > 0)
      all_pts(currents) = pts(currents);
    end
  end

  centers = centers(:, indexes);
  axes_length = axes_length(:, indexes);
  orient = orient(:, indexes);

  if (any(~isnan(axes_length(3, :))))
    axes_length = mymean(axes_length, 2);
  else
    mid_plane = folded_distribution(axes_length, sharpness, 2);
    z_coefs = get_struct('z-correlation');
    z_size = z_coefs.bkg + z_coefs.long_axis * mid_plane(1) + z_coefs.short_axis * mid_plane(2);

    axes_length = [mid_plane(1:2,1); z_size];
  end

  orient = align_orientations(orient);
  [axes_length, relative_z] = fit_relative_ellipse(all_pts, centers, axes_length, orient);

  centers = centers(:, indexes);
  orient = orient(1, indexes);
  [relative_z, centers, orient] = relative_position(pts, relative_z, centers, axes_length, orient);

  mymovie = align_embryo(mymovie, opts);

  if (abs(orient - mean(mymovie.dic.orientations)) > pi)
    orient = orient + pi;

    if (orient > 2*pi)
      orient = orient - 2*pi;
    end
  end

  mymovie.metadata.center_3d = centers;
  mymovie.metadata.axes_length_3d = axes_length;
  mymovie.metadata.orientation_3d = orient;
  mymovie.metadata.relative_z = relative_z;

  return;
end

function [relative_z, centers, orient] = relative_position(pts, z, centers, axes_length, orient)

  nframes = length(pts);
  relative_z = NaN(1, nframes);

  for i=1:nframes
    if (isnan(z(i)) & ~isempty(pts{i}))

      ell_coords = carth2elliptic(pts{i}(:, 1:2), centers(1:2, i), axes_length(1:2), orient(i), 'radial');
      dist = mymean(ell_coords(:,2));
      z_pos = axes_length(3)*sqrt(1 - dist.^2);

      if (imag(z_pos))

        sharpness = imnorm(pts{i}(:,5));
        val_thresh = graythresh(sharpness);

        goods = (sharpness >= val_thresh);
        goods = filter_goods(goods, 15);

        [tmp_c, tmp_a, tmp_o] = fit_ellipse(pts{i}(goods, 1:2));
        ell_coords = carth2elliptic(pts{i}(goods, 1:2), tmp_c, axes_length(1:2), tmp_o, 'radial');
        dist = mymean(ell_coords(:,2));
        z_pos = axes_length(3)*sqrt(1 - dist.^2);

        if (~imag(z_pos))
          bounding_box = [min(pts{i}(:,1:2)) max(pts{i}(:,1:2))];
          bb_center = [(bounding_box(1:2) + bounding_box(3:4))/2].';
          dist_thresh = [(bounding_box(3:4) - bounding_box(1:2))/4].';

          % Simple verification that the estimation did not fall completly off
          if (all(tmp_c >= bb_center - dist_thresh) && all(tmp_c <= bb_center + dist_thresh))

            relative_z(i) = z_pos;
            centers(1:2, i) = tmp_c;
            orient(i) = tmp_o;
          else
            display(['Weird estimation in frame ' num2str(i)])
          end
        end
      else
        relative_z(i) = z_pos;
      end
    else
      relative_z(i) = z(i);
    end
  end

  return;
end

function [axes_length, relative_z, all_coords] = fit_relative_ellipse(pts, centers, axes_length, orient)

  all_coords = NaN(0, size(pts{1}, 2)+1);
  neggs = length(pts);
  missing_frames = false(neggs, 1);

  for i=1:neggs
    if (~isempty(pts{i}))
      all_coords = [all_coords; [realign(pts{i}(:,1:2), [0;0], centers(1:2, i), orient(i)) pts{i}(:,3) - centers(3,i) pts{i}(:,4:5) ones(size(pts{i}(:,4)))*i]];
    else
      missing_frames(i) = true;
    end
  end

  pts = all_coords(:,1:3);
  npts = size(pts, 1);
  [frames, pindex] = unique(all_coords(:,6));
  z_pos = all_coords(pindex, 3);
  nframes = length(frames);
  indexes = all_coords(:,6);
  sharpness = imnorm(all_coords(:,5));
  thresh = graythresh(sharpness);

  unknowns = isnan(z_pos);

  w = 4;
  sharpness(sharpness < thresh) = sharpness(sharpness < thresh).^w;
  sharpness(sharpness >= thresh) = sharpness(sharpness >= thresh).^(1/w);

  all_errs = ones(nframes, 1);
  dist = NaN(nframes, 1);

  weights = NaN(size(sharpness));
  currents = false(npts, nframes);
  for i=1:nframes
    currents(:,i) = (indexes == frames(i));
    weights(currents(:,i)) = sharpness(currents(:,i)) / sum(sharpness(currents(:,i)));
  end

  z_coefs = get_struct('z-correlation');

  optims = optimset('Display', 'off', 'Algorithm', 'levenberg-marquardt', 'TolFun', 1e-6, 'TolX', 1e-4);

  p0 = sqrt(axes_length(1:2)*1.125);
  bests = lsqnonlin(@fit_z, p0, [], [], optims);

  bests = bests.^2;
  if (bests(1) < bests(2))
    [bests(1), bests(2)] = deal(bests(2), bests(1));
  end

  bests(3) = z_coefs.bkg + z_coefs.long_axis * bests(1) + z_coefs.short_axis * bests(2);
  axes_length = bests(1:3);

  [errs, rel_z] = fit_z(sqrt(axes_length));
  errs = errs(:);

  relative_z = NaN(1, neggs);
  relative_z(~missing_frames) = rel_z;

  return;

  function [z_err, tmp_z] = fit_z(ax)  

    ax = ax.^2;

    if (ax(1) < ax(2))
      [ax(1), ax(2)] = deal(ax(2), ax(1));
    end

    z_pred = z_coefs.bkg + z_coefs.long_axis * ax(1) + z_coefs.short_axis * ax(2);
    ax(3) = z_pred;

    z_err = all_errs;
    tmp_z = dist;
    z_std = dist;

    if (ax(3) < 0)
      z_err = 10;

      return;
    end

    z_target = sqrt(1 - (z_pos / ax(3)).^2);
    ell_coords = carth2elliptic(pts(:, 1:2), [0;0], ax(1:2), 0, 'radial');
    for j=1:nframes
      if (unknowns(j))
        z_dist = sum(ell_coords(currents(:,j), 2) .* weights(currents(:,j)));
      else
        z_dist = abs(sum(ell_coords(currents(:,j), 2) .* weights(currents(:,j))) - z_target(j));
      end

      if (~isfinite(z_dist) | z_dist > 1)
        z_err(j) = mean(sharpness(currents(:,j)));
      else
        if (unknowns(j))
          z_err(j) = sum(abs(ell_coords(currents(:, j), 2) - z_dist) .* weights(currents(:,j)));
          z_std(j) = sqrt(sum(((ell_coords(currents(:,j), 2) - z_dist).^2) .* weights(currents(:,j))));
          tmp_z(j) = ax(3)*sqrt(1 - z_dist^2) .* sign(z_dist);
        else
          z_err(j) = z_dist;
          z_std(j) = sqrt(sum((ell_coords(currents(:,j), 2) - z_target(j)).^2 .* weights(currents(:,j))) - z_dist^2);
          tmp_z(j) = z_pos(j);
        end
      end
    end

    if (all(isnan(tmp_z)))
      z_err = 10;

      return;
    end

    bads = (z_err > mean(z_err) + mymean(z_std));

    if (any(bads))
      tmp_z(bads) = NaN;

      all_bads = ismember(indexes, frames(bads));

      z_err(bads) = z_err(bads) .* mymean(sharpness(all_bads), 1, indexes(all_bads));
    end

    if (any(unknowns))
      z_err = mean(z_err) + 2*(sum(isnan(tmp_z))/nframes) + mymean(abs(tmp_z(unknowns)))/ax(3);
    else
      z_err = mean(z_err) + (sum(isnan(tmp_z))/nframes);
    end
    
    return;
  end
end

function [z_center, axes_length, pts] = fit_absolute_ellipse(pts, center, axes_length, orient)

  z_center = NaN;
  all_coords = NaN(0, size(pts{1}, 2));

  neggs = length(pts);
  for i=1:neggs
    all_coords = [all_coords; pts{i}];
  end
  all_coords = [realign(all_coords(:,1:2), [0;0], center, orient) all_coords(:,3:end)];

  npts = size(all_coords, 1);
  [planes, pindex] = unique(all_coords(:,4));
  z_pos = all_coords(pindex, 3);
  nplanes = length(planes);
  indexes = all_coords(:,4);
  sharpness = imnorm(all_coords(:,5));
  thresh = graythresh(sharpness);

  ell_coords = carth2elliptic(all_coords(:, 1:2), [0;0], axes_length, 0, 'radial');
  ell_coords(:,1) = ell_coords(:,1) - pi;
  [c, a, o] = fit_ellipse(all_coords(:, 3), ell_coords(:, 2) .* sign(ell_coords(:, 1)));
  z_center = c(1);

  w = 4;
  sharpness(sharpness < thresh) = sharpness(sharpness < thresh).^w;
  sharpness(sharpness >= thresh) = sharpness(sharpness >= thresh).^(1/w);

  all_errs = ones(nplanes, 1);
  dist = NaN(nplanes, 1);

  weights = NaN(size(sharpness));
  currents = false(npts, nplanes);
  for i=1:nplanes
    currents(:,i) = (indexes == planes(i));
    weights(currents(:,i)) = sharpness(currents(:,i)) / sum(sharpness(currents(:,i)));
  end

  z_coefs = get_struct('z-correlation');

  optims = optimset('Display', 'off', 'Algorithm', 'levenberg-marquardt', 'TolFun', 1e-6, 'TolX', 1e-4);

  p0 = [sqrt(axes_length(1:2)*1.125); z_center];
  bests = lsqnonlin(@fit_a, p0, [], [], optims);

  z_center = bests(3);

  bests = bests.^2;
  if (bests(1) < bests(2))
    [bests(1), bests(2)] = deal(bests(2), bests(1));
  end

  bests(3) = z_coefs.bkg + z_coefs.long_axis * bests(1) + z_coefs.short_axis * bests(2);
  axes_length = bests(1:3);

  return;

  function [z_err] = fit_a(ax)  

    c_z = z_pos - ax(3);
    ax = ax.^2;

    if (ax(1) < ax(2))
      [ax(1), ax(2)] = deal(ax(2), ax(1));
    end

    z_pred = z_coefs.bkg + z_coefs.long_axis * ax(1) + z_coefs.short_axis * ax(2);
    ax(3) = z_pred;

    z_err = all_errs;
    z_std = dist;

    if (ax(3) < 0)
      z_err = 10;

      return;
    end

    z_target = sqrt(1 - (c_z / ax(3)).^2);
    if (all(~isfinite(z_target)))
      z_err = 10;

      return;
    end

    ell_coords = carth2elliptic(all_coords(:, 1:2), [0;0], ax(1:2), 0, 'radial');
    for j=1:nplanes
      z_dist = abs(sum(ell_coords(currents(:,j), 2) .* weights(currents(:,j))) - z_target(j));

      if (~isfinite(z_dist) | z_dist > 1)
        z_err(j) = mean(sharpness(currents(:,j)));
        z_target(j) = NaN;
      else
        z_err(j) = z_dist;
        z_std(j) = sqrt(sum((ell_coords(currents(:,j), 2) - z_target(j)).^2 .* weights(currents(:,j))) - z_dist^2);
      end
    end

    bads = (z_err > mean(z_err) + mymean(z_std));

    if (any(bads))
      z_target(bads) = NaN;

      all_bads = ismember(indexes, planes(bads));

      z_err(bads) = z_err(bads) .* mymean(sharpness(all_bads), 1, indexes(all_bads));
    end

    z_err = mean(z_err) + (sum(isnan(z_target))/nplanes);
    
    return;
  end
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

%{
function plot_3d_ellipse(pts, axes, orient)
  
  hold on;

  if (nargin == 3)
    centers = pts;

    nslices = 7;
    z_pos = [-nslices:nslices] * (axes(3)) / (nslices+1) ;
    z_coef = sqrt(1 - (z_pos.^2 / axes(3).^2));
    z_pos = z_pos + centers(3);

    for i=1:length(z_pos)
      [x, y] = draw_ellipse(centers(1:2), axes(1:2)*z_coef(i), orient);
      plot3(x,y,z_pos(i)*ones(size(x)), 'r');
    end
  elseif (nargin == 2)
    z_pos = axes;
    if (iscell(pts))
      for i=1:length(pts)
        plot3(pts{i}(:,1),pts{i}(:,2),ones(size(pts{i},1),1)*z_pos(i), 'k');
      end
    else
      indxs = unique(pts(:,4));
      for i=1:length(indxs)
        currents = (pts(:,4) == indxs(i));
        plot3(pts(currents,1),pts(currents,2),ones(sum(currents), 1)*z_pos(i), 'k');
      end
    end
  else
    if (iscell(pts))
      for i=1:length(pts)
        plot3(pts{i}(:,1),pts{i}(:,2),pts{i}(:,3), 'k');
      end
    else
      plot3(pts(:,1),pts(:,2),pts(:,3), 'k');
    end
  end

  hold off;

  return;
end
%}
