function [best_coefs, shift] = calibrate_dic_correction(mymovie, trackings, opts)
%  Orig. Err ? std, Corr. Err ? std (% of Embryo radius)
%      0.0424    0.0312    0.0193    0.0174
%      a =
%          0.0424
%          -0.0440
%          -0.0176
%      b =
%
%          5.1836

  npts = 200;
  thetas = [0:2*pi/npts:2*pi]';
  thetas = thetas(1:end-1);

  if (ischar(mymovie))
    if (~isempty(findstr(mymovie, '*')))
      tmp = dir(mymovie); 
      mymovies = cell(1, length(tmp));

      intensity = [];
      corrections = [];
      ranges = [];
      for i=1:length(tmp)
        load(tmp(i).name);
        [tmp_intens, tmp_ranges, tmp_correction] = gather_corrections(mymovie, trackings, opts, thetas);
        intensity = cat(2, intensity, tmp_intens);
        ranges = cat(2, ranges, tmp_ranges);
        corrections = cat(2, corrections, tmp_correction);
      end
    else
      load(mymovie);
      [intensity, ranges, corrections] = gather_corrections(mymovie, trackings, opts, thetas);
    end
  else
    [intensity, ranges, corrections] = gather_corrections(mymovie, trackings, opts, thetas);
  end

  best_coefs = [];
  best_val = Inf;

  intensity = [intensity(:) ranges(:)];

  disp('can do better than this !')
  [shift] = fminbnd(@regress_dic, 0 , 2*pi);

  return;

  function curr_err = regress_dic(curr_shift)

    [~, indx] = min(abs(thetas-curr_shift))
    corrs = corrections([indx:end 1:indx-1], :, :);

    [coefs, stats] = robustfit(intensity, corrs(:));
    curr_err = stats.s;

    if (curr_err < best_val)
      best_coefs = coefs;
      best_val = curr_err;
    end

    return;
  end
end

function [intensity, ranges, corrections] = gather_corrections(mymovie, trackings, opts, thetas)

  fields = {};
  npts = length(thetas);

  %% Using the fluorescent channels as reference
  if (isempty(trackings))

    reference = mymovie.dic;
    target = mymovie.markers;

    if (isfield(mymovie, 'eggshell') & ~isempty(mymovie.eggshell))
      fields = [fields {'eggshell'}];
    end
    if (isfield(mymovie, 'cortex') & ~isempty(mymovie.cortex))
      fields = [fields {'cortex'}];
    end

  %% Using the trackings as reference
  else
    error('Not working yet!')
    
    reference = trackings.dic;
    markers = trackings.markers;
    
  end

  nframes = size_data(reference);

  centers = reference.centers;
  axes_length = reference.axes_length;
  orientations = reference.orientations;

  nfields = length(fields);

  egg_indx = -1;
  for f = 1:nfields
    if (strncmp(fields{f}, 'eggshell',8))
      egg_indx = f;
    end
  end

  dic_path = zeros(npts, nframes, nfields);
  fluo_path = zeros(npts, nframes, nfields);
  intensity = zeros(npts, nframes);
  maxval = double(intmax('uint16'));

  for nimg = 1:nframes
    img = double(load_data(mymovie.dic, nimg)) / maxval;

    dic_egg = reference.eggshell(nimg).carth;
    tmp_intens = bilinear(img, dic_egg(:,1), dic_egg(:,2));
    dic_egg = carth2elliptic(dic_egg, centers(:, nimg), axes_length(:, nimg), orientations(1, nimg));
    [~, intensity(:, nimg)] = interp_elliptic(dic_egg(:,1), tmp_intens, thetas);

    [~, dic_egg] = interp_elliptic(dic_egg, thetas);

    for f = 1:nfields
      tmp_path = target.(fields{f}).carth;
      tmp_path = carth2elliptic(tmp_path, centers(:, nimg), axes_length(:, nimg), orientations(1, nimg));
      [~, tmp_path] = interp_elliptic(tmp_path, thetas);

      fluo_path(:, nimg, f) = tmp_path;

      if (f == egg_indx)
        tmp_path = dic_egg;
      else
        tmp_path = reference.(fields{f}).carth;
        tmp_path = carth2elliptic(tmp_path, centers(:, nimg), axes_length(:, nimg), orientations(1, nimg));
        [~, tmp_path] = interp_elliptic(tmp_path, thetas);
      end

      dic_path(:, nimg, f) = tmp_path;
    end
  end

  intensity = repmat(intensity, [1 1 nfields]);
  mins = min(intensity);
  ranges = max(intensity) - mins; 
  ranges = repmat(ranges, [npts 1 1]);

  intensity = bsxfun(@minus, intensity, mins);
  intensity = intensity ./ ranges;

  corrections = dic_path - fluo_path;

  return;
end
