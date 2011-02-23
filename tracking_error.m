function trackings = tracking_error(trackings, opts)

  trackings = compute_reference(trackings, []);
  [trackings, paths] = recursive_error(trackings, opts);

  return;
end

function [trackings, paths] = recursive_error(trackings, opts)

  paths = cell(0);
  if (~isfield(trackings, 'child'))
    fields = fieldnames(trackings);

    for i=1:length(fields)
      field = fields{i};
      if (isstruct(trackings.(field)))
        trackings.(field) = recursive_error(trackings.(field), opts, hwait);
      end
    end

    return;
  end

  nchilds = length(trackings.child);
  if (isfield(trackings.child, 'shapes'))
    for i=1:nchilds
      [m,n,o] = size(trackings.child(i).shapes);
      paths(1:m,1:n,end+1:end+o) = trackings.child(i).shapes;
    end
  else
    for i=1:nchilds
      [trackings.child(i), tmp_paths] = recursive_error(trackings.child(i), opts);
      if (~isempty(tmp_paths))
        [m,n,o] = size(tmp_paths);
        paths(1:m,1:n,end+1:end+o) = tmp_paths;
      end
    end
  end

  trackings.errors = path_error(trackings.average, paths, trackings.reference, opts);

  return;
end

function mystruct = compute_reference(mystruct, orientations)
  
  if (isfield(mystruct, 'average'))
    [ntypes, nframes] = size(mystruct.average);

    if (isempty(mystruct.reference))
      mystruct.reference = get_struct('reference', 1);
    end
    for i=1:nframes
      [mystruct.reference.centers(:,i), mystruct.reference.axes_length(:,i), mystruct.reference.orientations(1,i)] = fit_ellipse(mystruct.average{1,i});
    end

    norients = length(orientations);
    orientations = [orientations mystruct.reference.orientations];
    new_orients = align_orientations(orientations);
    if (any(new_orients(1, 1:norients) ~= orientations(1, 1:norients)))
      orientations(1,norients+1:end) = orientations(1,norients+1:end) + pi;
    end
    orientations(1,orientations > 2*pi) = orientations(1,orientations > 2*pi) - 2*pi;

    mystruct.reference.orientations = orientations(1,norients+1:end);
  end

  if (isfield(mystruct, 'child'))
    for i = 1:length(mystruct.child)
      mystruct.child(i) = compute_reference(mystruct.child(i), orientations); 
    end
  end

  return;
end
