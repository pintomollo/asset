function mymovie = dp_data(mymovie, nimg, opts)

  if (isfield(mymovie, 'data'))
    if (isfield(mymovie.data, 'centers') && ~isempty(mymovie.data.centers))
      centers = mymovie.data.centers;
      axes_length = mymovie.data.axes_length;
      orientations = mymovie.data.orientations;
      eggshell = mymovie.data.eggshell;
      cortex = mymovie.data.cortex;

      if (~isfield(mymovie.data, 'update'))
        update = false(size(centers));
      else
        update = mymovie.data.update;
      end

      if (~isfield(mymovie.data, 'ruffles'))
        ruffles = get_struct('ruffles', size(eggshell));
      else
        ruffles = mymovie.data.ruffles;
      end
    else
      [nframes, imgsize] = size_data(mymovie.data);

      centers = NaN(2,nframes);
      axes_length = NaN(2,nframes);
      orientations = NaN(1,nframes);

      update = false(2,nframes);

      eggshell = get_struct('eggshell',[1,nframes]);
      cortex = get_struct('cortex',[1,nframes]);
      ruffles = get_struct('ruffles',[1, nframes]); 
    end

    mymovie.data.eggshell = eggshell;
    mymovie.data.cortex = cortex;
    mymovie.data.centers = centers;
    mymovie.data.axes_length = axes_length;
    mymovie.data.orientations = orientations;
    mymovie.data.update = update;
    mymovie.data.ruffles = ruffles;
  else
    error 'No suitable channels available for the data segmentation'; 

    return;
  end

  if (~isfield(mymovie.data, 'neighbors') | isempty(mymovie.data.neighbors) | empty_struct(mymovie.data.neighbors, 'axes_length'))
    mymovie = split_cells(mymovie, opts);
  end

  redo_nuclei = false;
  if (~isfield(mymovie.data, 'nuclei') | isempty(mymovie.data.nuclei) | empty_struct(mymovie.data.nuclei, 'carth'))
    if (all(isnan(mymovie.data.orientations)))
      redo_nuclei = true;
      mymovie.data.nuclei = get_struct('ruffles', [1,nframes]);
    else
      mymovie = detect_data_nuclei(mymovie, opts);
    end
  end

  neighbors = mymovie.data.neighbors;
  nuclei = mymovie.data.nuclei;
  parameters = opts.segmentation_parameters.data;

  img = [];

  if (length(cortex) < nimg | empty_struct(cortex(nimg), 'carth') | opts.recompute | (~strncmp(opts.do_ml, 'none', 4) & (strncmp(opts.ml_type, 'cortex', 6) | strncmp(opts.ml_type, 'all', 3))))

    update(2, nimg) = true;
    [centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), neighbors(nimg)] = detect_ellipse(neighbors(nimg), opts);

    if (isnan(orientations(1,nimg)))
      eggshell(nimg) = get_struct('eggshell');
      cortex(nimg) = get_struct('cortex');
      ruffles(nimg) = get_struct('ruffles');

      mymovie.data.centers = centers;
      mymovie.data.axes_length = axes_length;
      mymovie.data.orientations = orientations;
      mymovie.data.cortex = cortex;
      mymovie.data.eggshell = eggshell;
      mymovie.data.parameters = parameters;
      mymovie.data.update = update;
      mymovie.data.ruffles = ruffles;

      warning(['No embryo detected in frame ' num2str(nimg) ', skipping.']);
      return;
    end

    %img = imnorm(double(load_data(mymovie.data, nimg)));
    img = double(load_data(mymovie.data, nimg));
    noise_params = estimate_noise(img);
    vals = range(img(:));
    noise_thresh = min(round(vals/(8*noise_params(2)))+1, 25);
    signal_thresh = noise_params(2)*noise_thresh / vals;

    img = mask_neighbors(img, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), neighbors(nimg), opts);
      %img = mask_neighbors(img, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), neighbors(nimg), opts);

    if (opts.verbosity == 3)
      img_size = [330 450];
      figure;imshow(realign(imnorm(double(img)),img_size,centers(:,nimg),orientations(1,nimg)));
    end
  

    %polar_img = elliptic_coordinate(img, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), parameters.safety);

    %figure;imagesc(polar_img);

    %polar_img = imnorm(polar_img,[],[],'rows');

  
    %figure;imagesc(polar_img);

    %mask = roipoly(size(img,1),size(img,2), eggshell(nimg).carth(:,1), eggshell(nimg).carth(:,2));

    %noise = img;
    %noise(mask) = -1;
    %img = img .* mask;
    
    %minimg = mean(mean(noise(noise~=-1)));
    %img = imnorm(img, minimg, []);

    img = gaussian_mex(img, parameters.noise.gaussian);
    img = median_mex(img, parameters.noise.median);
    %img = imfilter(img,parameters.noise.gaussian,'symmetric');
    %img = medfilt2(img,parameters.noise.median);

    if (~isempty(nuclei(nimg).carth) & ~all(isnan(nuclei(nimg).carth)))
%      subplot(1,2,1);imagesc(img);
      img = img + GaussMask2D(nuclei(nimg).properties*0.75, size(img), nuclei(nimg).carth([2 1]), 0, 1)*prctile(img(:), 95);
%      subplot(1,2,2);imagesc(img);
%      drawnow;
    end

    img = imnorm(img, noise_params(1), []);

    %img = imnorm(img);

    polar_img = elliptic_coordinate(img, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), parameters.safety);
    
    %figure;imagesc(polar_img);
    %polar_img = imnorm(polar_img,[],[],'rows');

    polar_size = size(polar_img);

    %egg_path = carth2elliptic(eggshell(nimg).carth, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg));
    %egg_path = elliptic2pixels(egg_path, polar_size, axes_length(:,nimg), parameters.safety);

    %egg_path = adapt_path(polar_size, egg_path);

    %egg_path = adapt_path(size(polar_img), parameters.safety, ellpts);

    %parameters.cortex_weights.path = egg_path;

    if (opts.compute_probabilities)
      [cortex_path, emissions, transitions] = dynamic_programming(polar_img, parameters.cortex_params, parameters.scoring_func{2}, parameters.cortex_weights, opts);
      [cortex_path, emissions] = remove_polar_body(polar_img, cortex_path, parameters.cortex_params, parameters.scoring_func{2}, parameters.cortex_weights, opts, emissions);

      %'SET POWERS !!'
      [beta, gamma, probs] = find_temperatures(transitions, emissions, opts.temperatures);
      cortex(nimg).temperatures = [beta; gamma];
      %[entropy] = posterior_decoding(cortex_path, emissions, transitions, 1, 1);
    else

      parameters.cortex_params.nhood = 15;
      parameters.cortex_params.alpha = 0.1;
      parameters.cortex_params.beta = 0.85;
      parameters.cortex_params.gamma = 0.05;

      parameters.cortex_weights.alpha = 0.55;
      parameters.cortex_weights.beta = 0.025;
      %parameters.cortex_weights.gamma = 25*noise_params(2);
      parameters.cortex_weights.gamma = signal_thresh;

      %%%%%% Time-Lapse values
      %parameters.cortex_params.nhood = 11;
      %parameters.cortex_params.alpha = 0.4;
      %parameters.cortex_params.beta = 0.85;
      %parameters.cortex_params.gamma = 0.1;
      %parameters.cortex_weights.alpha = 0.55;
      %parameters.cortex_weights.beta = 0.05;
      %parameters.cortex_weights.gamma = 0.45;

      %print_all(parameters)
      %figure;imagesc(log(parameters.scoring_func{2}(polar_img, parameters.cortex_weights)));

      cortex_path = dynamic_programming(polar_img, parameters.cortex_params, parameters.scoring_func{2}, parameters.cortex_weights, opts);
      cortex_path = remove_polar_body(polar_img, cortex_path, parameters.cortex_params, parameters.scoring_func{2}, parameters.cortex_weights, opts);

%      figure;imshow(polar_img);
%      hold on;plot(cortex_path,[1:length(cortex_path)],'Color',[1 0.5 0]);

      %keyboard
    end

    ellpts = pixels2elliptic(cortex_path,size(polar_img), axes_length(:,nimg),parameters.safety);
    carths = elliptic2carth(ellpts,centers(:,nimg),axes_length(:,nimg),orientations(1,nimg));
    carths = carths(all(carths >= 1 & bsxfun(@le, carths , fliplr(size(img))), 2), :);

    if (opts.measure_performances)
      [estimation] = detect_ellipse(img, opts);
    end

    if (opts.verbosity == 3)

      figure;imagesc(parameters.scoring_func{2}(polar_img, parameters.cortex_weights))

      figure;imshow(polar_img);
      hold on;plot(cortex_path,[1:length(cortex_path)],'Color',[1 0.5 0]);
      %plot(egg_path,[1:length(cortex_path)],'g');

      figure;
      imshow(realign(img,img_size,centers(:,nimg),orientations(1,nimg)));
      hold on;
      myplot(realign(carths,img_size,centers(:,nimg),orientations(1,nimg)),'Color',[1 0.5 0]);
      %myplot(realign(eggshell(nimg).carth,img_size,centers(:,nimg),orientations(1,nimg)),'g');
      myplot(realign(draw_ellipse(centers(:,nimg), axes_length(:,nimg), orientations(1,nimg)),img_size,centers(:,nimg), orientations(1,nimg)),'m');

    end

    %cortex(nimg).raw = cortex_path;
    %cortex(nimg).elliptic = ellpts;
    cortex(nimg).carth = carths;

    if (opts.measure_performances)
      cortex(nimg).estim = estimation;
    end

    mymovie.data.centers = centers;
    mymovie.data.axes_length = axes_length;
    mymovie.data.orientations = orientations;
    mymovie.data.cortex = cortex;
  end

  if (length(eggshell) < nimg | empty_struct(eggshell(nimg), 'carth') | opts.recompute | (~strncmp(opts.do_ml, 'none', 4) & (strncmp(opts.ml_type, 'eggshell', 8) | strncmp(opts.ml_type, 'all', 3))))

    update(1, nimg) = true;

    parameters.eggshell_weights.alpha = 0.02;

    path = cortex(nimg).carth;
    indxs = convhull(path(:,1), path(:,2));
    path = path(indxs, :);

    [centers(:,nimg), axes_length(:,nimg), orientations(1,nimg)] = fit_ellipse(carths);
    egg_path = carth2elliptic(path, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), 'radial');

    max_rad = max(egg_path(:, 2)) * (1+parameters.eggshell_weights.alpha);
    axes_length(:, nimg) = axes_length(:, nimg) * max_rad;

    eggshell(nimg).carth = draw_ellipse(centers(:,nimg), axes_length(:,nimg), orientations(1,nimg));


    mymovie.data.centers = centers;
    mymovie.data.axes_length = axes_length;
    mymovie.data.orientations = orientations;
    mymovie.data.eggshell = eggshell;
  end

  if redo_nuclei
    opts.recompute = true;
    mymovie = detect_data_nuclei(mymovie, opts);
  end

  mymovie.data.parameters = parameters;
  mymovie.data.update = update;

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
  
  [center, axes_length, orientation] = deal(neighbors.centers(:, neighbors.index), ...
                                            neighbors.axes_length(:, neighbors.index), ...
                                            neighbors.orientations(:, neighbors.index));


  return;
end
