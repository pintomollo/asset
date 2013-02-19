function [datas, theta, ruffles] = gather_quantification(mymovie, opts)

%  opts = merge_structures(opts, get_struct('ASSET'));

  thresh = opts.quantification.pole_threshold;
  type = opts.quantification.kymograph_type;
  resolution = opts.quantification.resolution;

  nframes = size_data(mymovie.data);
  compute_ruffles = (nargout == 3);
  dist = NaN(nframes, 2);

  switch type
    case 'elliptic'
      npts = round(resolution);
      npts = npts + (mod(npts, 2) == 0);
      theta = [0:npts] / npts;
      theta = theta.' - 0.5;

      theta = theta * 2*pi;
      range = [-pi pi];
    case 'linear'
      npts = round(resolution);
      npts = npts + (mod((npts), 2) == 0);
      theta = [0:npts] / npts;
      theta = theta.' - 0.5;

      range = [-0.5 0.5];
    case 'projected'
      mymovie = reconstruct_egg(mymovie, opts);
  end

  indexes = [];
  tmp_pts = cell(nframes, 2 + compute_ruffles);

  for i=1:nframes
    warper = mymovie.data.warpers(i);

    if (isfield(mymovie.data.quantification, 'carth') & ~isempty(mymovie.data.quantification(i).carth))
      cortex = mymovie.data.quantification(i).carth;
    else
      cortex = mymovie.data.cortex(i).carth;
      if (opts.quantification.use_ruffles)
        cortex = insert_ruffles(cortex, mymovie.markers.ruffles(i).paths);
      end
    end

    if (compute_ruffles)
      ruffles_pts = ruffles_distance(cortex, mymovie.data.centers(:, i), mymovie.data.axes_length(:, i), mymovie.data.orientations(i));
    end

    switch type
      case 'elliptic'
        cortex = carth2normalized(cortex, warper, opts);
        pts = carth2elliptic(cortex, warper.reference.centers, warper.reference.axes_length, warper.reference.orientations);
        pts = pts(:, 1);

        if (any(pts(2:end) < pts(1:end-1)))
          goods = true(size(pts));
          prev = -Inf;
          for p=1:length(pts)
            if (pts(p) < prev)
              goods(p) = false;
            else
              prev = pts(p);
            end
          end

          indexes = [1:length(pts)];
          indexes = indexes(goods);

          pts = pts(goods);
        else
          indexes = [];
        end
      case 'linear'
        cortex = carth2normalized(cortex, warper, opts);
        pts = carth2linear(cortex);
      case 'direct'
        cortex = realign(cortex, [0;0], mymovie.data.centers(:, i), mymovie.data.orientations(i));
        cortex = cortex * opts.pixel_size;

        [pts, indexes, dist(i, :)] = realign_poles(cortex, thresh, opts);
      case 'normalized'
        cortex = carth2normalized(cortex, warper, opts);

        [pts, indexes, dist(i, :)] = realign_poles(cortex, thresh, opts);
      case 'projected'
        if (isnan(mymovie.metadata.relative_z(i)))
          pts = [];
        else
          cortex = cortex * opts.pixel_size;
          cortex = project2midplane(cortex, mymovie.metadata.center_3d(:, i), mymovie.metadata.axes_length_3d, mymovie.metadata.orientation_3d(i), mymovie.metadata.relative_z(i));
          cortex = realign(cortex, [0;0], mymovie.metadata.center_3d(1:2, i), mymovie.metadata.orientation_3d(i));

          [pts, indexes, dist(i, :)] = realign_poles(cortex, thresh, opts);
        end
    end
    tmp_pts{i,1} = pts;
    tmp_pts{i,2} = indexes;

    if (compute_ruffles)
      tmp_pts{i,3} = ruffles_pts;
    end
  end

  switch type
    case {'elliptic', 'linear'}
      dist(:, 1) = range(1);
      dist(:, 2) = range(2);

      full_indexes = theta;
    case {'direct', 'normalized', 'projected'}
      boundaries = max(abs(dist));

      boundaries(1) = -boundaries(1);
      full_indexes = [fliplr([0:-resolution:boundaries(1)]) resolution:resolution:boundaries(2)];
      npts = length(full_indexes);

  end

  if (npts > 1e4)
    warning('Too many elements necessary for the quantification');
    datas = [];
    ruffles = [];
    theta = [];

    return;
  end

  datas = NaN(nframes, npts);
  if (compute_ruffles)
    ruffles = NaN(nframes, npts);
  end

  for i=1:nframes

    pts = tmp_pts{i,1};
    indexes = tmp_pts{i,2};

    if (isempty(pts))
      continue;
    end

    if (isempty(indexes))
      new_intens = interp_elliptic(pts, mymovie.data.quantification(i).cortex, full_indexes, dist(i,:));
    else
      new_intens = interp_elliptic(pts, mymovie.data.quantification(i).cortex(indexes), full_indexes, dist(i,:));
    end
    new_intens(full_indexes < dist(i,1) | full_indexes > dist(i,2),:) = NaN;
    datas(i, :) = new_intens(:, 2);

    if (compute_ruffles)
      ruffles_pts = tmp_pts{i, 3};

      if (isempty(indexes))
        new_ruffles = interp_elliptic(pts, ruffles_pts, full_indexes, dist(i,:));
      else
        new_ruffles = interp_elliptic(pts, ruffles_pts(indexes), full_indexes, dist(i,:));
      end
      new_ruffles(full_indexes < dist(i,1) | full_indexes > dist(i,2), :) = NaN;
      ruffles(i, :) = new_ruffles(:, 2);
    end
  end

  theta = full_indexes;

  goods = ~all(isnan(datas), 2);

  if (any(~goods))
    first = find(goods, 1, 'first');
    last = find(goods, 1, 'last');
    goods([1:first last:end]) = true;

    [X,Y] = meshgrid(theta, [1:size(datas, 1)].');
    datas(~goods, :) = interp2(X(goods, :), Y(goods, :), datas(goods, :), X(~goods, :), Y(~goods, :), 'linear');
  end

  if (nargout == 3)
    tmp = theta;
    theta = ruffles;
    ruffles = tmp;
  end

  return;
end

function [pts, indexes, dists] = realign_poles(cortex, thresh, opts)

  [pts, total_dist] = carth2linear(cortex, opts);
  
  goods = (cortex(:,1) > 0);
  max_pos = max(cortex(goods, 1));
  goods = (cortex(:,1) > (1-thresh)*max_pos);

  max_r = min(abs(cortex(goods,2)));
  post_indx = find(goods & abs(cortex(:,2)) == max_r, 1);

  goods = (cortex(:,1) < 0);
  max_pos = min(cortex(goods, 1));
  goods = (cortex(:,1) < (1-thresh)*max_pos);

  max_r = min(abs(cortex(goods,2)));
  ant_indx = find(goods & abs(cortex(:,2)) == max_r, 1);

  pts = pts * total_dist;
  pts = pts - pts(post_indx);
  if (ant_indx > post_indx)
    pts = [pts(ant_indx:end, :) - total_dist; pts(1:ant_indx-1, :)];
  else
    pts = [pts(ant_indx:end, :); pts(1:ant_indx-1, :) + total_dist];
  end

  dists = [pts(1), total_dist + pts(1)];
  indexes = [ant_indx:length(pts) 1:ant_indx-1];

  return;
end

function dist = ruffles_distance(pts, center, axes_length, orient)

  conv = pts(convhull(pts(:,1),pts(:,2)),:);
  if (~ispolycw(conv(:,1),conv(:,2)))
    conv = conv([end:-1:1],:);
  end

  ell_pts = carth2elliptic(pts, center, axes_length, orient, 'radial');
  ell_conv = carth2elliptic(conv, center, axes_length, orient, 'radial');

  [junk, indx] = sort(ell_conv(:, 1));
  ell_conv = ell_conv(indx, :);
  [ell_pos, indx] = sort(ell_pts(:,1));
  [junk, back_indx] = sort(indx);

  tmp_conv = ell_conv([end 1:end 1],:);
  tmp_conv(1,1) = tmp_conv(1,1) - 2*pi;
  tmp_conv(end,1) = tmp_conv(end,1) + 2*pi;
  tmp_conv = interp1q(tmp_conv(:,1), tmp_conv(:,2), ell_pos);
  tmp_conv = tmp_conv(back_indx);

  ell_conv = [ell_pts(:,1) tmp_conv];
  dist = ell_conv(:,2) - ell_pts(:,2);

  return;
end
