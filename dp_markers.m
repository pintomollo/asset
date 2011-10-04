function mymovie = dp_markers(mymovie, nimg, opts)

  if (isfield(mymovie, 'markers') && isfield(mymovie.markers, 'centers') && ~isempty(mymovie.markers.centers))
    centers = mymovie.markers.centers;
    axes_length = mymovie.markers.axes_length;
    orientations = mymovie.markers.orientations;
    eggshell = mymovie.markers.eggshell;
    cortex = mymovie.markers.cortex;
    if (~isfield(mymovie.dic, 'neighbors'))
      neighbors = get_struct('reference', size(eggshell));
    else
      neighbors = mymovie.dic.neighbors;
    end

    if (~isfield(mymovie.markers, 'update'))
      update = false(size(centers));
    else
      update = mymovie.markers.update;
    end

  elseif (isfield(mymovie, 'cortex') && isfield(mymovie, 'eggshell'))
    [nframes, imgsize] = size_data(mymovie.eggshell);

    centers = zeros(2,nframes);
    axes_length = zeros(2,nframes);
    orientations = zeros(1,nframes);

    update = false(2,nframes);

    eggshell = get_struct('eggshell',[1,nframes]);
    cortex = get_struct('cortex',[1,nframes]);
    neighbors = get_struct('reference',[1, nframes]);
  else
    error 'No suitable channels available for the markers segmentation'; 

    return;
  end

  parameters = opts.segmentation_parameters.markers;

  img = [];
  if (~isempty(mymovie.eggshell) & ( length(eggshell) < nimg | isempty(eggshell(nimg).carth) | opts.recompute | (~strncmp(opts.do_ml,'none',4) & strncmp(opts.ml_type, 'eggshell', 8))))

    update(1, nimg) = true;

    img = imnorm(double(load_data(mymovie.eggshell,nimg)));

    if (opts.verbosity == 3)
      figure;
      imshow(realign(img,[388 591],centers(:,nimg),orientations(1,nimg)));
    end

    img = gaussian_mex(img, parameters.noise.gaussian);
    img = median_mex(img, parameters.noise.median);
    %img = imfilter(img, parameters.noise.gaussian,'symmetric');
    %img = medfilt2(img,parameters.noise.median);
    img = imnorm(img);
    img = 1-img;

    if (opts.measure_performances)
      [centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), neighbors(nimg), estimation] = detect_ellipse(img, false, opts);
    else
      [centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), neighbors(nimg)] = detect_ellipse(img, false, opts);
    end

    if (opts.verbosity == 3)
      old_center = centers(:,nimg);
      old_axes = axes_length(:,nimg);
      old_orient = orientations(1,nimg);
    end

    img = mask_neighbors(img, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), neighbors(nimg), opts);

    polar_img = elliptic_coordinate(img, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), parameters.safety);
    polar_img = imnorm(polar_img,[],[],'rows');

    if (opts.compute_probabilities)
      [egg_path, emissions, transitions] = dynamic_programming(polar_img, parameters.eggshell_params, parameters.scoring_func{1}, parameters.eggshell_weights, opts);

      %'SET POWERS !!'
      %[entropy] = posterior_decoding(egg_path, emissions, transitions, 1, 1);
      [beta, gamma, probs] = find_temperatures(transitions, emissions, opts.temperatures);
      eggshell(nimg).temperatures = [beta; gamma];
    else
      egg_path = dynamic_programming(polar_img, parameters.eggshell_params, parameters.scoring_func{1}, parameters.eggshell_weights, opts);
    end

    ellpts = pixels2elliptic(egg_path,size(polar_img), axes_length(:,nimg),parameters.safety);
    carths = elliptic2carth(ellpts,centers(:,nimg),axes_length(:,nimg),orientations(1,nimg));
    carths = carths(all(carths > 0 & bsxfun(@le, carths , fliplr(size(img))), 2), :);

    [centers(:,nimg), axes_length(:,nimg), orientations(1,nimg)] = fit_ellipse(carths);

    %eggshell(nimg).raw = egg_path;
    %eggshell(nimg).elliptic = ellpts;
    eggshell(nimg).carth = carths;

    if (opts.verbosity == 3)
      %hax = gca;
      %image(img, 'Parent', hax);
      figure;
      imshow(realign(img,[388 591],centers(:,nimg),orientations(1,nimg)));
      hold on;
      myplot(realign(carths,[388 591],centers(:,nimg),orientations(1,nimg)),'g');
      myplot(realign(draw_ellipse(old_center, old_axes, old_orient),[388 591],old_center, old_orient),'m');
      myplot(realign(draw_ellipse(centers(:,nimg), axes_length(:,nimg), orientations(1,nimg)),[388 591],centers(:,nimg), orientations(1,nimg)),'c');

      figure;
      imshow(polar_img);
      hold on;
      plot(egg_path, [1:length(egg_path)], 'g');
      %ell = elliptic2pixels([0 1; 2*pi 1], size(polar_img), parameters.safety);
      %plot(ell(:,2), [1 length(egg_path)], 'm');
      plot(ones(1,2) * size(polar_img,2) * 5/ 6, [1 length(egg_path)], 'm');
    end

    if (opts.measure_performances)
      eggshell(nimg).estim = estimation;
    end

    mymovie.markers.centers = centers;
    mymovie.markers.axes_length = axes_length;
    mymovie.markers.orientations = orientations;
    mymovie.markers.eggshell = eggshell;
    mymovie.markers.neighbors = neighbors;
  end

  
  if (~isempty(mymovie.cortex) & (length(cortex) < nimg | isempty(cortex(nimg).carth) | opts.recompute | (~strncmp(opts.do_ml,'none',4) & strncmp(opts.ml_type, 'cortex', 6))))

    update(2, nimg) = true;
    img = imnorm(double(load_data(mymovie.cortex,nimg)));
    img = mask_neighbors(img, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), neighbors(nimg), opts);

    if (opts.verbosity == 3)
      figure;imshow(realign(img,[388 591],centers(:,nimg),orientations(1,nimg)));
    end

    mask = roipoly(size(img,1),size(img,2), eggshell(nimg).carth(:,1), eggshell(nimg).carth(:,2));

    noise = img;
    noise(mask) = -1;;
    img = img .* mask;
    
    minimg = mean(mean(noise(noise~=-1)));
    img = imnorm(img, minimg, []);

    img = gaussian_mex(img, parameters.noise.gaussian);
    img = median_mex(img, parameters.noise.median);
    %img = imfilter(img,parameters.noise.gaussian,'symmetric');
    %img = medfilt2(img,parameters.noise.median);
    img = imnorm(img);

    polar_img = elliptic_coordinate(img, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), parameters.safety);
    polar_size = size(polar_img);

    egg_path = carth2elliptic(eggshell(nimg).carth, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg));
    egg_path = elliptic2pixels(egg_path, polar_size, axes_length(:,nimg), parameters.safety);

    egg_path = adapt_path(polar_size, egg_path);

    %egg_path = adapt_path(size(polar_img), parameters.safety, ellpts);

    parameters.cortex_weights.path = egg_path;

    if (opts.compute_probabilities)
      [cortex_path, emissions, transitions] = dynamic_programming(polar_img, parameters.cortex_params, parameters.scoring_func{2}, parameters.cortex_weights, opts);
      [cortex_path, emissions] = remove_polar_body(polar_img, cortex_path, parameters.cortex_params, parameters.scoring_func{2}, parameters.cortex_weights, opts, emissions);

      %'SET POWERS !!'
      [beta, gamma, probs] = find_temperatures(transitions, emissions, opts.temperatures);
      cortex(nimg).temperatures = [beta; gamma];
      %[entropy] = posterior_decoding(cortex_path, emissions, transitions, 1, 1);
    else
      cortex_path = dynamic_programming(polar_img, parameters.cortex_params, parameters.scoring_func{2}, parameters.cortex_weights, opts);
      cortex_path = remove_polar_body(polar_img, cortex_path, parameters.cortex_params, parameters.scoring_func{2}, parameters.cortex_weights, opts);
    end

    ellpts = pixels2elliptic(cortex_path,size(polar_img), axes_length(:,nimg),parameters.safety);
    carths = elliptic2carth(ellpts,centers(:,nimg),axes_length(:,nimg),orientations(1,nimg));
    carths = carths(all(carths > 0 & bsxfun(@le, carths , fliplr(size(img))), 2), :);

    if (opts.measure_performances)
      [estimation] = detect_ellipse(img, true, opts);
    end

    if (opts.verbosity == 3)
      figure;imshow(polar_img);
      hold on;plot(cortex_path,[1:length(cortex_path)],'Color',[1 0.5 0]);
      plot(egg_path,[1:length(cortex_path)],'g');

      figure;
      imshow(realign(img,[388 591],centers(:,nimg),orientations(1,nimg)));
      hold on;
      myplot(realign(carths,[388 591],centers(:,nimg),orientations(1,nimg)),'Color',[1 0.5 0]);
      myplot(realign(eggshell(nimg).carth,[388 591],centers(:,nimg),orientations(1,nimg)),'g');
      myplot(realign(draw_ellipse(centers(:,nimg), axes_length(:,nimg), orientations(1,nimg)),[388 591],centers(:,nimg), orientations(1,nimg)),'m');

    end

    %cortex(nimg).raw = cortex_path;
    %cortex(nimg).elliptic = ellpts;
    cortex(nimg).carth = carths;

    if (opts.measure_performances)
      cortex(nimg).estim = estimation;
    end

    mymovie.markers.cortex = cortex;
  end

  mymovie.markers.parameters = parameters;
  mymovie.markers.update = update;

  return;
end

function [center, axes_length, orientation, neighbors, estim] = detect_ellipse(img, estim_only, opts)

  center = NaN(2, 1);
  axes_length = NaN(2, 1);
  orientation = NaN;
  estim = [];
  neighbors = get_struct('reference');

  neighbors.centers = center;
  neighbors.axes_length = axes_length;
  neighbors.orientations = orientation;

  imgsize = size(img);
  npixels = max(imgsize);

  size250 = round(npixels/250);
  size200 = round(npixels/200);

  img = gaussian_mex(img, size250);
  img = median_mex(img, size200, 3);
  img = imnorm(img);

  thresh = graythresh(img);
  img = (img > thresh*(max(img(:))) );

  [ellipse, estim] =  split_cells(img, estim_only, opts);

  if (estim_only)
    center = estim;
    axes_length = [];
    orientation = [];
    estim = [];

    return;
  end

  if (numel(ellipse) == 0)
    return;
  elseif (size(ellipse, 1) > 1)
    imgsize = size(img);

    %dist = sum(bsxfun(@minus, ellipse(:, [1 2]), imgsize([2 1])/2).^2, 2);

    dist = [bsxfun(@minus, ellipse(:, [1 2]), imgsize([2 1])).^2 ellipse(:, [1 2]).^2];
    dist = min(dist, [], 2);

    [~, indx] = max(dist);
    indxs = false(size(dist));
    indxs(indx(1)) = true;

    [neighbors.centers, neighbors.axes_length, neighbors.orientations] = deal(ellipse(~indxs, 1:2).', ellipse(~indxs, 3:4).', ellipse(~indxs, 5));

    ellipse = ellipse(indxs,:);
  end

  [center, axes_length, orientation] = deal(ellipse(1:2).', ellipse(3:4).', ellipse(5));

  return;
end
