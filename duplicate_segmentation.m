function mymovie = duplicate_segmentation(mymovie, type, opts, nframe)

  switch (opts.segmentation_type)
    case 'dic'
      [nframes imgsize ] = size_data(mymovie.dic);
      orig_type = 'dic';
    case {'markers', 'all'}
      if (isfield(mymovie, 'eggshell') & ~isempty(mymovie.eggshell))
        [nframes imgsize ] = size_data(mymovie.eggshell);
      else
        [nframes imgsize ] = size_data(mymovie.cortex);
      end
      orig_type = 'markers';
    otherwise
      error 'None of the expected field are present in ''mymovie''';
  end

  if (isfield(mymovie, type) && isfield(mymovie.(type), 'centers'))
    centers = mymovie.(type).centers;
    axes_lengths = mymovie.(type).axes_length;
    orientations = mymovie.(type).orientations;
    eggshell = mymovie.(type).eggshell;
    cortex = mymovie.(type).cortex;
  else
    centers = zeros(2,nframes);
    axes_lengths = zeros(2,nframes);
    orientations = zeros(1,nframes);
    eggshell = get_struct('eggshell',[1,nframes]);
    cortex = get_struct('cortex',[1,nframes]);
  end

  %mymovie.(type) = mymovie.dic;
  if (isfield(mymovie.(orig_type), 'inverted') & ~isempty(mymovie.(orig_type).inverted))
    mymovie.(type).inverted = mymovie.(orig_type).inverted;
  else
    mymovie.(orig_type).inverted = false;
    mymovie.(type).inverted = false;
  end

  if (isfield(mymovie.(orig_type), 'cytokinesis'))
    mymovie.(type).cytokinesis = mymovie.(orig_type).cytokinesis;
  end

  if (nargin == 4)
    frames = nframe;
  else
    frames = 1:nframes;
  end

  if (isfield(mymovie.(orig_type), 'neighbors'))
    neighbors = mymovie.(orig_type).neighbors;
  else
    neighbors = get_struct('reference', [1 nframes]);
  end

  if (strncmp(orig_type, 'dic', 3))
    params = opts.segmentation_parameters.correction;
    maxval = double(intmax('uint16'));

    for i=1:length(frames)
      nimg = frames(i);

      img = double(load_data(mymovie.dic,nimg)) / maxval;

      egg = mymovie.dic.eggshell(nimg).carth;
      cor = mymovie.dic.cortex(nimg).carth;
      intens = bilinear_mex(img, egg);
      mins = min(intens);
      range = max(intens) - mins;
      intens = (intens - mins) / range;

      center = mymovie.dic.centers(:,nimg);
      axes_length = mymovie.dic.axes_length(:,nimg);
      orientation = mymovie.dic.orientations(:,nimg);

      egg = carth2elliptic(egg, center, axes_length, orientation) + params.safety;
      cor = carth2elliptic(cor, center, axes_length, orientation);

      [egg, indx] = align(egg);
      cor = align(cor);
      intens = align(intens, indx);

      theta = egg(:,1);
      shift = abs(theta - params.shift);
      indx = find(shift == min(shift), 1);
      intens = intens([end-indx+2:end 1:end-indx+1]);

      correction = (params.bkg + intens*params.factor + range*params.range);
      new_egg = egg(:,2) - correction;

      correction = correction([end 1:end 1]);
      theta = [theta(end)-2*pi; theta(1:end); theta(1)+2*pi];
      new_cortex = cor(:,2) - interp1q(theta, correction, cor(:,1));

      new_egg = elliptic2carth(egg(:,1), new_egg, center, axes_length, orientation);
      new_cortex = elliptic2carth(cor(:,1), new_cortex, center, axes_length, orientation);

      eggshell(nimg).carth = new_egg;
      cortex(nimg).carth = new_cortex;
      [centers(:,nimg), axes_lengths(:,nimg), orientations(1,nimg)] = fit_ellipse(new_egg);

      neighbors(nimg).centers = bsxfun(@plus, neighbors(nimg).centers, centers(:,nimg) - center);
      neighbors(nimg).orientations = neighbors(nimg).orientations + (orientations(:,nimg) - orientation);
    end
    orientations = align_orientations(orientations, mean(mymovie.dic.orientations));
  else
    for i=1:length(frames)
      nimg = frames(i);
      eggshell(nimg).carth = mymovie.(orig_type).eggshell(nimg).carth;
      cortex(nimg).carth = mymovie.(orig_type).cortex(nimg).carth;
      centers(:, nimg) = mymovie.(orig_type).centers(:, nimg);
      axes_lengths(:, nimg) = mymovie.(orig_type).axes_length(:, nimg);
      orientations(:, nimg) = mymovie.(orig_type).orientations(:, nimg);
    end
  end

  mymovie.(type).eggshell = eggshell;
  mymovie.(type).cortex = cortex;
  mymovie.(type).centers = centers;
  mymovie.(type).axes_length = axes_lengths;
  mymovie.(type).orientations = orientations;
  mymovie.(type).neighbors = neighbors;

  return;
end

function [ellpts, indx] = align(ellpts, indx)

  if (nargin < 2)
    indx = find(ellpts(:,1)==min(ellpts(:,1)));
  end

  if (all(ellpts(1,:)) == ellpts(end,:))
    ellpts = [ellpts(indx:end-1,:); ellpts(1:indx-1,:); ellpts(indx,:)];
  else
    ellpts = [ellpts(indx:end,:); ellpts(1:indx-1,:)];
  end

  return;
end
