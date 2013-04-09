function mymovie = correct_jitter(mymovie, opts)

  type = opts.segmentation_type;
  if (strncmp(type, 'all', 3))
    opts.segmentation_type = 'dic';
    mymovie = correct_jitter(mymovie, opts)
    opts.segmentation_type = 'markers';
    mymovie = correct_jitter(mymovie, opts)

    return;
  end

  mymovie.(type) = align_orientations(mymovie.(type));

  orientations = mymovie.(type).orientations;
  axes_length = mymovie.(type).axes_length;

  valids = all(isfinite([orientations; axes_length]), 1);
  orientations = orientations(valids);
  axes_length = axes_length(:, valids);

  eccentricity = sqrt(1 - (axes_length(2,:).^2) ./ (axes_length(1,:).^2));
  %confidence = eccentricity;
  confidence = (erf(5*(eccentricity - 0.5)) + 1) / 2;

  aim = median(orientations);
  dorient = orientations - aim;

%  figure;
%  hist(dorient, 20);

  dorient = dorient .* confidence;
  orientations = aim + dorient;

  mymovie.(type).orientations(valids) = orientations;

  if (isfield(mymovie, 'metadata') && isfield(mymovie.metadata, 'orientation_3d'))

    orientations = mymovie.metadata.orientation_3d;
    axes_length = mymovie.metadata.axes_length_3d;

    valids = isfinite(orientations);
    orientations = orientations(valids);

    eccentricity = sqrt(1 - (axes_length(2).^2) ./ (axes_length(1).^2));
    confidence = (erf(5*(eccentricity - 0.5)) + 1) / 2;

    aim = median(orientations);
    dorient = orientations - aim;

    dorient = dorient .* confidence;
    orientations = aim + dorient;

    mymovie.metadata.orientation_3d(valids) = orientations;
  end

  return;
end
