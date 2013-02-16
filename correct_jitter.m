function mymovie = correct_jitter(mymovie, opts)

  type = opts.segmentation_type;
  if (strncmp(type, 'all', 3))
    opts.segmentation_type = 'dic';
    mymovie = correct_jitter(mymovie, opts)
    opts.segmentation_type = 'markers';
    mymovie = correct_jitter(mymovie, opts)

    return;
  end

  orientations = mymovie.(type).orientations;
  axes_length = mymovie.(type).axes_length;

  eccentricity = sqrt(1 - (axes_length(2,:).^2) ./ (axes_length(1,:).^2));
  %confidence = eccentricity;
  confidence = (erf(5*(eccentricity - 0.5)) + 1) / 2;

  aim = median(orientations);
  dorient = orientations - aim;

%  figure;
%  hist(dorient, 20);

  dorient = dorient .* confidence;
  orientations = aim + dorient;

  mymovie.(type).orientations = orientations;

  if (isfield(mymovie, 'metadata') && isfield(mymovie.metadata, 'orientation_3d'))

    orientations = mymovie.metadata.orientation_3d;
    axes_length = mymovie.metadata.axes_length_3d;

    eccentricity = sqrt(1 - (axes_length(2).^2) ./ (axes_length(1).^2));
    confidence = (erf(5*(eccentricity - 0.5)) + 1) / 2;

    aim = median(orientations);
    dorient = orientations - aim;

    dorient = dorient .* confidence;
    orientations = aim + dorient;

    mymovie.metadata.orientation_3d = orientations;
  end

  return;
end
