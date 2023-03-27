function [trackings] = mean_trackings(trackings, opts, hwait)

  if (nargin < 3)
    if (nargin == 1)
      opts = get_struct('ASSET');
    end
    if (opts.verbosity > 1)
      hwait = waitbar(0, 'Computing tracking means', 'Name', 'ASSET');
    else
      hwait = [];
    end
    nsteps = 1;
    nindx = 0;
  elseif (ishandle(hwait))
    nsteps = get(hwait, 'UserData');
    nindx = nsteps(2);
    nsteps = nsteps(1);
  else
    nsteps = 1;
    nindx = 0;
  end

  means = {[]};

  if (~isfield(trackings, 'child'))
    fields = fieldnames(trackings);

    for i=1:length(fields)
      field = fields{i};
      if (isstruct(trackings.(field)))
        trackings.(field) = mean_trackings(trackings.(field), opts, hwait);
      end
    end

    return;
  end

  nchilds = length(trackings.child);
  nsteps = nsteps * nchilds;

  for i=1:nchilds
    if (isfield(trackings.child(i), 'shapes'))
      if (opts.follow_periphery)
        tmp_paths = path_periphery(trackings.child(i).shapes);
      else
        tmp_paths = trackings.child(i).shapes;
      end

      [m,n] = size(tmp_paths);
      means(1:m,1:n,i) = tmp_paths;
    else
      if (opts.verbosity > 1)
        set(hwait, 'UserData', [nsteps (nindx + (i-1)/nsteps)])
      end
      trackings.child(i) = mean_trackings(trackings.child(i), opts, hwait);

      [m,n] = size(trackings.child(i).average);
      means(1:m,1:n,i) = trackings.child(i).average;
    end

    if (opts.verbosity > 1)
      waitbar(nindx + (i / nsteps), hwait);
    end
  end

  if (nchilds > 1)
    trackings.average = mean_paths(means, opts);
  else
    trackings.average = means;
  end
  
  if (nargin < 3 & opts.verbosity > 1)
    close(hwait);
  end

  return;
end

