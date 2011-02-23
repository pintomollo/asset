function mymovie = dp_dic(mymovie, parameters, nimg, opts)

  if (isfield(mymovie, 'dic'))
    if (isfield(mymovie.dic, 'centers'))
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
    else
      [nframes, imgsize] = size_data(mymovie.dic);

      centers = zeros(2,nframes);
      axes_length = zeros(2,nframes);
      orientations = zeros(1,nframes);

      update = false(2,nframes);

      eggshell = get_struct('eggshell',[1,nframes]);
      cortex = get_struct('cortex',[1,nframes]);
    end
  else
    error 'No DIC available for the DIC segmentation'; 

    return;
  end

  img = [];
  global rescale_size;
  rescale_size = [388 591]*1;

  if (length(eggshell) < nimg | isempty(eggshell(nimg).carth) | opts.recompute | (~strncmp(opts.do_ml, 'none', 4) & strncmp(opts.ml_type, 'eggshell', 8)))
    update(1,nimg) = true;

    img=load_data(mymovie.dic,nimg);

    if (opts.measure_performances)
      [centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), estimation] = detect_ellipse(img, opts.verbosity, false);
    else
      [centers(:,nimg), axes_length(:,nimg), orientations(1,nimg)] = detect_ellipse(img, opts.verbosity, false);
    end

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

    [egg_path, egg_shift] = detect_eggshell(polar_img, outer_egg, axes_length(:,nimg), parameters.eggshell_weights.eta);

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

    ell_egg = pixels2elliptic(egg_path,size(polar_img),parameters.safety);
    cart_egg = elliptic2carth(ell_egg,centers(:,nimg),axes_length(:,nimg),orientations(1,nimg));
  
    [centers(:,nimg), axes_length(:,nimg), orientations(1,nimg)] = fit_ellipse(cart_egg);
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

  
  if (length(cortex) < nimg | isempty(cortex(nimg).carth) | opts.recompute | (~strncmp(opts.do_ml, 'none', 4) & strncmp(opts.ml_type, 'cortex', 6)))
    update(2,nimg) = true;

    if (isempty(img))
      img = load_data(mymovie.dic, nimg);
    end

    polar_img = elliptic_coordinate(img, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), parameters.safety);
    polar_size = size(polar_img);

    egg_path = carth2elliptic(eggshell(nimg).carth, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg));
    egg_path = elliptic2pixels(egg_path, polar_size, parameters.safety);
    egg_path = adapt_path(polar_size, egg_path);

    [inner_egg, outer_egg] = detect_eggshell(polar_size, egg_path, axes_length(:,nimg), parameters.eggshell_weights.eta, eggshell(nimg).thickness);
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
      estimation = detect_ellipse(eggless, opts.verbosity, true);

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

    ellpts = pixels2elliptic(cortex_path,size(polar_img),parameters.safety);
    carths = elliptic2carth(ellpts,centers(:,nimg),axes_length(:,nimg),orientations(1,nimg));

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
      ell = elliptic2pixels([0 1; 2*pi 1], polar_size, parameters.safety);
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
  end

  mymovie.dic.parameters = parameters;
  mymovie.dic.update = update;

  return;
end

function [center, axes_length, orientation, estim] = detect_ellipse(img, verbosity, estim_only)

  global rescale_size;

  if (nargin == 1)
    estim_only = false;
    verbosity = 0;
  end

  npixels = max(size(img));
  size10 = round(npixels/10);
  size75 = round(npixels/75);
  size150 = round(npixels/150);
  size100 = round(prod(size(img)) / 100);

  img = imadm_mex(img);

  thresh = graythresh(img);
  edg1 = (img>thresh*0.5);

  if (verbosity == 3)
    old_edg1 = edg1;
  end

  edg1 = increase_canevas(edg1, size10);

  edg1 = imdilate(edg1, strel('disk', size150));
  edg1 = imfill(edg1, 'holes');
  edg1 = bwareaopen(edg1, size100);
  edg1 = imdilate(edg1, strel('disk', size75));
  edg1 = imfill(edg1, 'holes');
  edg1 = imerode(edg1, strel('disk', size75 + size150));

  edg1 = reduce_canevas(edg1, size10);

  if(~any(any(edg1)))
    display('Error : No embryo detected !!');
  end

  estim = bwboundaries(edg1, 8, 'noholes');
  estim = estim{1};
  estim = estim(:,[2 1]);

  if (estim_only)
    center = estim;
    axes_length = [];
    orientation = [];
    estim = [];

    return
  end

  [center, axes_length, orientation] = estimate_ellipse(edg1);
  if (verbosity == 3)
    figure; imshow(realign(img,rescale_size,center,orientation));
    figure; imshow(realign(old_edg1,rescale_size,center,orientation));
    figure; imshow(realign(edg1,rescale_size,center,orientation));
    hold on; myplot(realign(draw_ellipse(center, axes_length, orientation),rescale_size,center, orientation),'m');
  end

  return
end

%***********************************************************
function ellsize = elliptic_size(mymovie,safety)

  axes_length = ceil(mymovie.axes_length * safety);
  perifs = ceil(ellipse_circum(axes_length))-1;
  ellsize = [max(perifs) max(axes_length(1,:))];
  
end

%***********************************************************
function resized = increase_canevas(img, new_size)

  [h,w] = size(img);
  
  resized =  zeros(h+(2*new_size), w+(2*new_size));

  if (islogical(img))
    resized = logical(resized);
  end

  resized((new_size+1):(new_size+h),(new_size+1):(new_size+w)) = img;
  
  return;
end

%***********************************************************
function resized = reduce_canevas(img, new_size)

  [h,w] = size(img);
  
  resized =  img((new_size+1):(h-new_size), (new_size+1):(w-new_size));

  return;
end
