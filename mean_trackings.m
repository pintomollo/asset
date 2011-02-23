function [trackings] = mean_trackings(trackings, opts, hwait)
  
  if (nargin < 3)
    if (opts.verbosity > 0)
      hwait = waitbar(0, 'Computing tracking means', 'Name', 'RECOS');
    else
      hwait = [];
    end
    nsteps = 1;
    nindx = 0;
  elseif (ishandle(hwait))
    nsteps = get(hwait, 'UserData');
    nindx = nsteps(2);
    nsteps = nsteps(1);
  end

  mean_splines = get_struct('spline',[0 0 0]);

  if (isfield(trackings, 'child'))
    nchilds = length(trackings.child);
  else
    nchilds = 0;
  end

  if (isfield(trackings, 'files'))
    nfiles = length(trackings.files);
  else
    nfiles = 0;
  end

  nsteps = nsteps * (nchilds + nfiles);

  for i=1:nchilds
    set(hwait, 'UserData', [nsteps (nindx + (i-1)/nsteps)])
    trackings.child(i) = mean_trackings(trackings.child(i), opts, hwait);
    [m,n] = size(trackings.child(i).mean);
    mean_splines(1:m,1:n,i) = trackings.child(i).mean;

    waitbar(nindx + (i / nsteps), hwait);
  end

  for i=1:nfiles
    if (isempty(trackings.files(i).splines) || opts.recompute)
      trackings.files(i).splines = shapes2splines(trackings.files(i).shapes);

      if (opts.follow_periphery)
        trackings.files(i).splines = spline_periphery(trackings.files(i).splines);
      end
    end
    [m,n] = size(trackings.files(i).splines);
    mean_splines(1:m,1:n,i+nchilds) = trackings.files(i).splines;

    waitbar(nindx + (i + nchilds) / nsteps, hwait);
  end

  if (nchilds + nfiles > 1)
    trackings.mean = mean_spline(mean_splines, opts);
  else
    trackings.mean = mean_splines;
  end
  
  if (nargin < 3)
    close(hwait);
  end

  return;
end

