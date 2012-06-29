function [datas, flows] = cortical_flow(mymovie, opts)

  resolution = 0.5;
  thresh = 1/10;

  mymovie = duplicate_segmentation(mymovie, 'dic', opts);
  nframes = size_data(mymovie.dic);

  dpos = [-20:1:0];
  nstrip = length(dpos);

  tmp_pts = cell(nframes, 2);
  dist = NaN(nframes, 2);

  for i=1:nframes
    nimg = i;

    img = double(load_data(mymovie.dic, nimg));
    %cortex = insert_ruffles(mymovie.dic.cortex(nimg).carth, mymovie.dic.ruffles(nimg).paths);
    cortex = mymovie.dic.cortex(nimg).carth;

    values = perpendicular_sampling(img, cortex, dpos, opts);

    warper = mymovie.markers.warpers(nimg);
    cortex = carth2normalized(cortex, warper, opts);
    [pts, total_dist] = carth2linear(cortex);
    npts = length(pts);

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

    dist(i,1) = pts(1);
    dist(i,2) = total_dist + pts(1);

    tmp_pts{i,1} = pts;
    tmp_pts{i,2} = values([ant_indx:npts 1:ant_indx-1], :);
  end

  boundaries = max(abs(dist));
  if (any(boundaries > 1e4))
    keyboard
    error('Too many elements');
  end

  boundaries(1) = -boundaries(1);
  full_indexes = [fliplr([0:-resolution:boundaries(1)]) resolution:resolution:boundaries(2)];
  npts = length(full_indexes);
  datas = NaN(nframes, npts, nstrip);
  flows = datas;

  for i=1:nframes
    pts = tmp_pts{i,1}; % * opts.pixel_size;
    [junk, new_intens] = interp_elliptic(pts, tmp_pts{i, 2}, full_indexes, dist(i,:));
    new_intens(full_indexes < dist(i,1) | full_indexes > dist(i,2),:) = NaN;
    %datas([1:nstrip] + (i-1)*nstrip, :) = new_intens.';
    datas(i, :, :) = reshape(new_intens, [1 npts, nstrip]);
  end

  params = get_struct('smoothness_parameters');
  params.nhood = 5;
  params.alpha = 1;
  params.beta = 0.25;
  params.gamma = 0.25;

  func = @(x)(ones(size(x)));
  for i=1:nstrip
    [junk, flows(:,:,i)] = dynamic_programming(datas(:,:,i), params, func, [], opts, true);
  end

  return;
end
