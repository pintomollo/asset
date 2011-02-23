function trackings = tracking_error(trackings, opts)

  trackings = compute_reference(trackings, []);
  [trackings, splines] = recursive_error(trackings, opts);

  return;
end

function [trackings, splines] = recursive_error(trackings, opts)

  splines = get_struct('spline',[0 0 0]);
  if (isfield(trackings, 'child'))
    nchilds = length(trackings.child);
    for i=1:nchilds
      [trackings.child(i), tmp_splines] = recursive_error(trackings.child(i), opts);
      if (~isempty(tmp_splines))
        [m,n,o] = size(tmp_splines);
        splines(1:m,1:n,end+1:end+o) = tmp_splines;
      end
    end
  else
    nchilds = 0;
  end

  for i=1:length(trackings.files)
    [m,n,o] = size(trackings.files(i).splines);
    splines(1:m,1:n,end+1:end+o) = trackings.files(i).splines;
  end

  trackings.errors = spline_error(trackings, splines, opts.nbins);

  return;
end

function mystruct = compute_reference(mystruct, orientations)
  
  if (isfield(mystruct, 'mean'))
    [ntypes, nframes] = size(mystruct.mean);

    if (isempty(mystruct.reference))
      mystruct.reference = get_struct('reference', 1);
    end
    for i=1:nframes
      [mystruct.reference.center(:,i), mystruct.reference.axes_length(:,i), mystruct.reference.orientation(1,i)] = fit_ellipse(mystruct.mean(1,i));
    end


    norients = length(orientations);
    orientations = [orientations mystruct.reference.orientation];
    new_orients = align_orientations(orientations);
    if (any(new_orients(1, 1:norients) ~= orientations(1, 1:norients)))
      orientations(1,norients+1:end) = orientations(1,norients+1:end) + pi;
    end
    orientations(1,orientations > 2*pi) = orientations(1,orientations > 2*pi) - 2*pi;

    mystruct.reference.orientation = orientations(1,norients+1:end);

    mystruct.elliptic = get_struct('spline', [ntypes, nframes]);
    for i=1:ntypes
      for j=1:nframes
        mystruct.elliptic(i,j) = carth2elliptic(mystruct.mean(i,j), mystruct.reference.center(:,j), mystruct.reference.axes_length(:,j), mystruct.reference.orientation(1,j));
      end
    end
  end

  if (isfield(mystruct, 'child'))
    for i = 1:length(mystruct.child)
      mystruct.child(i) = compute_reference(mystruct.child(i), orientations); 
    end
  end

  return;
end
