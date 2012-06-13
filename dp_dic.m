function mymovie = dp_dic(mymovie, nimg, opts)

  if (isfield(mymovie, 'dic'))
    if (isfield(mymovie.dic, 'centers') & ~isempty(mymovie.dic.centers))
      centers = mymovie.dic.centers;
      axes_length = mymovie.dic.axes_length;
      orientations = mymovie.dic.orientations;
      eggshell = mymovie.dic.eggshell;
      cortex = mymovie.dic.cortex;

      if (~isfield(mymovie.dic, 'update'))
        update = false(size(centers));
      else
        update = mymovie.dic.update;
      end

      if (~isfield(mymovie.dic, 'ruffles'))
        ruffles = get_struct('ruffles', size(eggshell));
      else
        ruffles = mymovie.dic.ruffles;
      end
    else
      [nframes, imgsize] = size_data(mymovie.dic);

      centers = NaN(2,nframes);
      axes_length = NaN(2,nframes);
      orientations = NaN(1,nframes);

      update = false(2,nframes);

      eggshell = get_struct('eggshell',[1,nframes]);
      cortex = get_struct('cortex',[1,nframes]);
      ruffles = get_struct('ruffles',[1, nframes]);

      if (~isfield(mymovie.dic, 'neighbors') | isempty(mymovie.dic.neighbors))
        mymovie = split_cells(mymovie, opts);
      end
    end

    mymovie.dic.eggshell = eggshell;
    mymovie.dic.cortex = cortex;
    mymovie.dic.centers = centers;
    mymovie.dic.axes_length = axes_length;
    mymovie.dic.orientations = orientations;
    mymovie.dic.update = update;
    mymovie.dic.ruffles = ruffles;

  else
    error 'No DIC available for the DIC segmentation'; 

    return;
  end

  %keyboard
 %nimg

  parameters = opts.segmentation_parameters.dic;
      mymovie.dic.parameters = parameters;

  img = [];
  global rescale_size;
  rescale_size = [388 591]*1;

  if (length(eggshell) < nimg | isempty(eggshell(nimg).carth) | opts.recompute | (~strncmp(opts.do_ml, 'none', 4) & (strncmp(opts.ml_type, 'eggshell', 8) | strncmp(opts.ml_type, 'all', 3))))
    update(1,nimg) = true;

    img = imnorm(double(load_data(mymovie.dic,nimg)));

    %figure;imshow(img);hold on;
    %title(num2str(nimg));

    if (opts.measure_performances)
      [centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), neighbors(nimg), estimation] = detect_ellipse(mymovie.dic.neighbors(nimg), img, opts);
    else
    %beep;keyboard
      [centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), neighbors(nimg)] = detect_ellipse(mymovie.dic.neighbors(nimg), opts);
    end
    if (isnan(orientations(1,nimg)))
      mymovie.dic.eggshell(nimg).warped = [];

      mymovie.dic.eggshell(nimg) = get_struct('eggshell');
      mymovie.dic.centers = centers;
      mymovie.dic.axes_length = axes_length;
      mymovie.dic.orientations = orientations;
      mymovie.dic.cortex(nimg) = get_struct('cortex');
      %mymovie.dic.neighbors(nimg) = get_struct('reference');

      warning(['No embryo detected in frame ' num2str(nimg) ', skipping.']);
      return;
    end

    img = mask_neighbors(img, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), neighbors(nimg), opts);
    %orientations(1,nimg) = orientations(1,nimg) + pi;

    polar_img = elliptic_coordinate(img, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), parameters.safety);

    if (opts.compute_probabilities)
      [outer_egg, emissions, transitions] = dynamic_programming(polar_img, parameters.eggshell_params, parameters.scoring_func{1}, parameters.eggshell_weights, opts);

      [beta, gamma, probs] = find_temperatures(transitions, emissions, opts.temperatures);

      if (opts.verbosity == 3)
        figure;imshow(imnorm(probs));colormap(jet)
      end

      probs = carthesian_coordinate(probs,centers(:,nimg),axes_length(:,nimg),orientations(1,nimg),parameters.safety,size(img),[false true]);
      probs(isnan(probs)) = 0;
      eggshell(nimg).temperatures = [beta; gamma];
      eggshell(nimg).probabilities = sparse(probs);

      if (opts.verbosity == 3)
        figure;imagesc(realign(probs,rescale_size,centers(:,nimg),orientations(1,nimg)));
      end
    else
      outer_egg = dynamic_programming(polar_img, parameters.eggshell_params, parameters.scoring_func{1}, parameters.eggshell_weights, opts);
    end

    [egg_path, egg_shift] = detect_eggshell(polar_img, outer_egg, axes_length(:,nimg), parameters.safety, parameters.eggshell_weights.eta);
    %keyboard
    if (~(egg_shift > 0))
      mymovie.dic.eggshell(nimg).warped = [];

      mymovie.dic.eggshell(nimg) = get_struct('eggshell');
      mymovie.dic.centers(:, nimg) = NaN;
      mymovie.dic.axes_length(:, nimg) = NaN;
      mymovie.dic.orientations(1, nimg) = NaN;
      %mymovie.dic.neighbors(nimg) = neighbors;
      mymovie.dic.cortex(nimg) = get_struct('cortex');

      warning(['No embryo detected in frame ' num2str(nimg) ', skipping.']);
      return;
    end


    if (opts.verbosity == 3)
      figure; imshow(polar_img)
      hold on; plot(outer_egg, [1:length(outer_egg)], 'g');
      %ell = elliptic2pixels([0 1; 2*pi 1], size(polar_img), parameters.safety);
      %plot(ell(:,2), [1 length(egg_path)], 'm');
      plot(ones(1,2) * size(polar_img,2) * 5/ 6, [1 length(egg_path)], 'm');

      old_center = centers(:,nimg);
      old_axes = axes_length(:,nimg);
      old_orient = orientations(1,nimg);

      figure;imshow(imadm_mex(polar_img));
    end

    ell_egg = pixels2elliptic(egg_path,size(polar_img),axes_length(:,nimg),parameters.safety);
    cart_egg = elliptic2carth(ell_egg,centers(:,nimg),axes_length(:,nimg),orientations(1,nimg));

    cart_egg = cart_egg(all(cart_egg >= 1 & bsxfun(@le, cart_egg , fliplr(size(img))), 2), :);

    [centers(:,nimg), axes_length(:,nimg), orientations(1,nimg)] = fit_ellipse(cart_egg);
    %plot(cart_egg(:,1), cart_egg(:,2), 'g');
    %draw_ellipse(centers(:, nimg), axes_length(:,nimg), orientations(1,nimg), 'r');

    %orientations(1,nimg) = orientations(1,nimg) + pi

    %eggshell(nimg).raw = [inner_egg egg_path outer_egg];
    %eggshell(nimg).elliptic = ell_egg;
    eggshell(nimg).carth = cart_egg;
    eggshell(nimg).thickness = egg_shift;

    if (opts.measure_performances)
      eggshell(nimg).estim = estimation;
    end
    %if (opts.compute_probabilities)
    %  eggshell(nimg).entropy = entropy;
    %end
    if (opts.verbosity == 3)
      figure;
      imshow(realign(img,rescale_size,centers(:,nimg),orientations(1,nimg)));
      hold on;
      myplot(realign(cart_egg,rescale_size,centers(:,nimg),orientations(1,nimg)),'g');
      myplot(realign(draw_ellipse(old_center, old_axes, old_orient),rescale_size,old_center, old_orient),'m');
      myplot(realign(draw_ellipse(centers(:,nimg), axes_length(:,nimg), orientations(1,nimg)),rescale_size,centers(:,nimg), orientations(1,nimg)),'c');
    end

    mymovie.dic.eggshell = eggshell;
    mymovie.dic.centers = centers;
    mymovie.dic.axes_length = axes_length;
    mymovie.dic.orientations = orientations;
  end

  if (length(cortex) < nimg | isempty(cortex(nimg).carth) | opts.recompute | (~strncmp(opts.do_ml, 'none', 4) & (strncmp(opts.ml_type, 'cortex', 6) | strncmp(opts.ml_type, 'all', 3))))
    update(2,nimg) = true;

    if (isempty(img))
      img = imnorm(double(load_data(mymovie.dic, nimg)));
      img = mask_neighbors(img, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), neighbors(nimg), opts);
    end

    polar_img = elliptic_coordinate(img, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), parameters.safety);
    polar_size = size(polar_img);

    egg_path = carth2elliptic(eggshell(nimg).carth, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg));

    if (isempty(egg_path))
      return;
    end

    egg_path = elliptic2pixels(egg_path, polar_size, axes_length(:,nimg), parameters.safety);
    egg_path = adapt_path(polar_size, egg_path);

    [inner_egg, outer_egg] = detect_eggshell(polar_size, egg_path, axes_length(:,nimg), parameters.safety, parameters.eggshell_weights.eta, eggshell(nimg).thickness);
    %inner_egg = adapt_path(size(polar_img), parameters.safety, eggshell(nimg).raw(:,1));
    %egg_path = adapt_path(size(polar_img), parameters.safety, eggshell(nimg).raw(:,2));
    %outer_egg = adapt_path(size(polar_img), parameters.safety, eggshell(nimg).raw(:,3));

    parameters.cortex_weights.path = [inner_egg egg_path outer_egg];
    %parameters.cortex_weights.egg = egg_path;
    %parameters.cortex_weights.inner = inner_egg;
    %parameters.cortex_weights.outer = outer_egg;

    polar_cortex = erase_egg(polar_img, outer_egg, inner_egg);

    if (opts.measure_performances)
      eggless = carthesian_coordinate(polar_cortex,centers(:,nimg),axes_length(:,nimg),orientations(1,nimg),parameters.safety,size(img),[false true]);
      estimation = detect_ellipse(eggless, opts);

      if (opts.verbosity == 3)
        figure;imshow(realign(eggless,rescale_size,centers(:,nimg),orientations(1,nimg)));
        myplot(realign(estimation,rescale_size,centers(:,nimg),orientations(1,nimg)));
      end
    end

    if (opts.compute_probabilities)
      [cortex_path, emissions, transitions] = dynamic_programming(polar_cortex, parameters.cortex_params, parameters.scoring_func{2}, parameters.cortex_weights, opts);
      %[cortex_path, emissions] = remove_polar_body(polar_img, cortex_path, parameters.cortex_params, parameters.scoring_func{2}, parameters.cortex_weights, opts, emissions);

      %'SET POWERS !!'
      [beta, gamma, probs] = find_temperatures(transitions, emissions, opts.temperatures);
      cortex(nimg).temperatures = [beta; gamma];
      probs = carthesian_coordinate(probs,centers(:,nimg),axes_length(:,nimg),orientations(1,nimg),parameters.safety,size(img),[false true]);
      probs(isnan(probs)) = 0;
      cortex(nimg).probabilities = sparse(probs);
      %[entropy] = posterior_decoding(cortex_path, emissions, transitions, 2, 2)
    else
      cortex_path = dynamic_programming(polar_cortex, parameters.cortex_params, parameters.scoring_func{2}, parameters.cortex_weights, opts);
      %cortex_path = remove_polar_body(polar_img, cortex_path, parameters.cortex_params, parameters.scoring_func{2}, parameters.cortex_weights, opts);
    end

    ellpts = pixels2elliptic(cortex_path,size(polar_img), axes_length(:,nimg),parameters.safety);
    carths = elliptic2carth(ellpts,centers(:,nimg),axes_length(:,nimg),orientations(1,nimg));

    carths = carths(all(carths >= 1 & bsxfun(@le, carths , fliplr(size(img))), 2), :);

    %cortex(nimg).raw = cortex_path;
    %cortex(nimg).elliptic = ellpts;
    cortex(nimg).carth = carths;

    if (opts.measure_performances)
      cortex(nimg).estim = estimation;
    end
    %if (opts.compute_probabilities)
    %  cortex(nimg).entropy = entropy;
    %end
    if (opts.verbosity == 3)
      figure;imshow(polar_img);
      hold on;plot(inner_egg,[1:length(inner_egg)]);
      plot(outer_egg,[1:length(outer_egg)]);

      figure;imshow(polar_cortex);
      hold on;plot(cortex_path,[1:length(cortex_path)],'Color',[1 0.5 0]);
      plot(egg_path,[1:length(cortex_path)],'g');
      ell = elliptic2pixels([0 1; 2*pi 1], polar_size, axes_length(:,nimg), parameters.safety);
      plot(ell(:,2), [1 length(cortex_path)], 'm');
      %plot(carths(:,1),carths(:,2),'r');

      figure;
      imshow(realign(img,rescale_size,centers(:,nimg),orientations(1,nimg)));
      hold on;
      myplot(realign(carths,rescale_size,centers(:,nimg),orientations(1,nimg)),'Color',[1 0.5 0]);
      myplot(realign(eggshell(nimg).carth,rescale_size,centers(:,nimg),orientations(1,nimg)),'g');
      myplot(realign(draw_ellipse(centers(:,nimg), axes_length(:,nimg), orientations(1,nimg)),rescale_size,centers(:,nimg), orientations(1,nimg)),'m');
    end

    mymovie.dic.cortex = cortex;
    mymovie.dic.ruffles(nimg) = get_struct('ruffles');
  end

  mymovie.dic.parameters = parameters;
  mymovie.dic.update = update;

  return;
end

function [center, axes_length, orientation, neighbors, estim] = detect_ellipse(neighbors, img, opts)

  if (nargin == 2)
    if (isstruct(neighbors))
      opts = img;
      img = [];
      estim_only = false;
    else
      tmp = neighbors;
      opts = img;
      img = tmp;
      neighbors = [];
      estim_only = true;
    end
  else
    estim_only = false;
  end

  center = NaN(2, 1);
  axes_length = NaN(2, 1);
  orientation = NaN;
  estim = [];

  if (~isempty(img))
    [ellipse, estim] =  split_cells(img, true, opts);
  end

  if (estim_only)
    center = estim;
    axes_length = [];
    orientation = [];
    estim = [];

    return;
  end

  if (all(isnan(neighbors.centers(:))) | neighbors.index <= 0)
    return;
  end

  %elseif (size(neighbours.centers, 2) > 1)
    %imgsize = size(img);

    %dist = sum(bsxfun(@minus, ellipse(:, [1 2]), imgsize([2 1])/2).^2, 2);

    %dist = [bsxfun(@minus, ellipse(:, [1 2]), imgsize([2 1])).^2 ellipse(:, [1 2]).^2];
    %dist = min(dist, [], 2);

    %[~, indx] = max(dist);
    %indxs = false(size(dist));
    %indxs(indx(1)) = true;

    %[neighbors.centers, neighbors.axes_length, neighbors.orientations] = deal(ellipse(~indxs, 1:2).', ellipse(~indxs, 3:4).', ellipse(~indxs, 5));

  %  ellipse = ellipse(indxs,:);
  %end

  [center, axes_length, orientation] = deal(neighbors.centers(:, neighbors.index), ...
                                            neighbors.axes_length(:, neighbors.index), ...
                                            neighbors.orientations(:, neighbors.index));

  return;
end
