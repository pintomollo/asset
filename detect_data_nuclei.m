function mymovie = detect_data_nuclei(mymovie, opts)

  [nframes, imgsize] = size_data(mymovie.data);
  parameters = opts.segmentation_parameters.data;

  if (isfield(mymovie.data, 'nuclei') & ~isempty(mymovie.data.nuclei))
    nuclei = mymovie.data.nuclei;
  else
    nuclei = get_struct('ruffles', [1,nframes]);
  end

  if (~isfield(opts, 'nuclei_tracking'))
    opts.nuclei_tracking = get_struct('spot_tracking');
    opts = load_parameters(opts, 'nuclei_tracking');
  end

  if (opts.recompute | empty_struct(nuclei, 'carth'))

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

    radius_thresh = 0.5 / opts.pixel_size;
    %surf_thresh = 0.5;
    surf_thresh = 0.1;

    for i=1:nframes
      nimg = i;
      %nimg = 290;
      %nimg = randi(nframes, 1);

      nuclei(nimg).cluster = [];
      nuclei(nimg).carth = [];
      nuclei(nimg).properties = [];

      if (~isnan(mymovie.data.orientations(1,nimg)))

        img = double(load_data(mymovie.data, nimg));

        if (opts.verbosity == 3)
          orig_img = img;
        end

        noise_params = estimate_noise(img);
        vals = range(img(:));

        inner_thresh = min(round(vals / (15*noise_params(2))) + 1, 35);
        noise_thresh = min(round(vals / (15*noise_params(2))) + 1, 15);

        img = gaussian_mex(img, parameters.noise.gaussian);
        img = median_mex(img, parameters.noise.median);
        threshed = (img > noise_params(1) + noise_thresh*noise_params(2));
        bw = imfill(threshed, 'holes');

        cc = bwconncomp(bw);
        if (cc.NumObjects > 1)
          sizes = NaN(1, cc.NumObjects);
          for c = 1:cc.NumObjects
            sizes(c) = numel(cc.PixelIdxList{c});
          end
          [v, indx] = max(sizes);
        else
          indx = 1;
        end
        bw = false(size(img));
        bw(cc.PixelIdxList{indx}) = true;
        pix = cc.PixelIdxList{indx};
        [x,y] = ind2sub(size(bw), pix);
        indx = convhull(x,y);

        threshed = (img > noise_params(1) + inner_thresh*noise_params(2));

        bw = roipoly(bw, y(indx), x(indx));
        bw = imerode(bw, strel('disk', 5));

        if (opts.verbosity == 3)
          bkg_mask = bw;
        end

        bw = bw & ~threshed;

        cc = bwconncomp(bw);
        dist = bwdist(~bw);

        %props = regionprops(bw, 'Area', 'Eccentricity', 'Centroid', 'EquivDiameter', 'Solidity');

        if (opts.verbosity == 3)
          img_size = [330 450]
          figure;imagesc(realign(img, img_size, mymovie.data.centers(:,nimg), mymovie.data.orientations(1,nimg)+pi));

          figure;imagesc(realign(double(bkg_mask + threshed), img_size, mymovie.data.centers(:,nimg), mymovie.data.orientations(1,nimg)+pi));

          figure;imagesc(realign(bw, img_size, mymovie.data.centers(:,nimg), mymovie.data.orientations(1,nimg)+pi));
          hold on;
        end

        for j = 1:cc.NumObjects
          vals = dist(cc.PixelIdxList{j});
          [R, indx] = max(vals);

          if (R > radius_thresh & (2*pi*R^2 / numel(vals)) > surf_thresh)

            indx = cc.PixelIdxList{j}(indx(1));
            [y, x] = ind2sub(size(bw), indx);
        
            if (opts.verbosity == 3)
              tmp_pts = realign([x y], img_size, mymovie.data.centers(:,nimg), mymovie.data.orientations(1,nimg)+pi);
              rectangle('Position', [tmp_pts(1)-R tmp_pts(2)-R 2*R([1 1])], 'Curvature', [1 1], 'EdgeColor', 'w');
            end

            nuclei(nimg).carth = [nuclei(nimg).carth; [x y]];
            nuclei(nimg).properties(end+1, 1) = R;
          end
        end
      end

      text = sprintf('Progress: %d/%d', i, nframes);
      cpb.setValue(i);
      cpb.setText(text);
    end

    cpb.stop();
  end

%  figure;imagesc(load_data(mymovie.data, 100));
%  hold on
%  for i=1:nframes
%    if (~isempty(nuclei(i).carth))
%      scatter(nuclei(i).carth(:,1), nuclei(i).carth(:,2));
%    end
%  end
  

  if (opts.recompute | empty_struct(nuclei, 'cluster'))
    opts.nuclei_tracking.pixel_size = opts.pixel_size;
    nuclei = track_spots(nuclei, opts.nuclei_tracking);
  end

  path_length = cell(nframes, 1);
  max_length = 0;
  for i=1:nframes
    nimg = i;

    links = nuclei(nimg).cluster;
    path_length{nimg} = zeros(size(nuclei(nimg).carth, 1), 1);
    for j=1:size(links, 1)
      path_length{nimg}(links(j,1)) = path_length{links(j,3)}(links(j,2)) + nimg - links(j,3);
    end
    if (any(path_length{nimg} > max_length))
      max_length = max(path_length{nimg});
    end
  end
  all_pos = NaN(nframes,3);
  for i=nframes:-1:1
    nimg = i;

    links = nuclei(nimg).cluster
    for j=1:size(links, 1)
      path_length{links(j,3)}(links(j,2)) = path_length{nimg}(links(j,1));
    end

    goods = (path_length{nimg} == max_length);
    if (any(goods))
      all_pos(nimg, :) = [nuclei(nimg).carth(goods,:) nuclei(nimg).properties(goods)];
    end
  end

  if (max_length < opts.nuclei_tracking.min_path_length)
    goods = false(nframes, 1);
    warning('No nucleus detected !');
  else
    indxs = [1:nframes].';
    goods = ~isnan(all_pos(:,1));

    smooth_level = opts.nuclei_tracking.smoothing_complexity - 1;
    if (smooth_level >= 0)
      tmp = emdc(indxs(goods), all_pos(goods,1));
      all_pos(goods,1) = sum(tmp(max(end-smooth_level,1):end,:), 1);
      tmp = emdc(indxs(goods), all_pos(goods,2));
      all_pos(goods,2) = sum(tmp(max(end-smooth_level,1):end,:), 1);
      tmp = emdc(indxs(goods), all_pos(goods,3));
      all_pos(goods,3) = sum(tmp(max(end-smooth_level,1):end,:), 1);
    end

    if (opts.nuclei_tracking.interpolate)
      first = find(goods, 1, 'first');
      last = find(goods, 1, 'last');

      all_pos(first:last, :) = interp1q(indxs(goods), all_pos(goods, :), indxs(first:last));
    end
  end

  prev_good = 0;
  for i=1:nframes
    if (goods(i))
      nuclei(i).carth = all_pos(i,1:2);
      nuclei(i).properties = all_pos(i,3);
      if (prev_good > 0)
        nuclei(i).cluster = [1 1 prev_good];
      end

      prev_good = i;
    else
      nuclei(i).carth = NaN(0,2);
      nuclei(i).properties = [];
      nuclei(i).cluster = NaN(0,3);
    end
  end

  mymovie.data.nuclei = nuclei;

  return;
end
