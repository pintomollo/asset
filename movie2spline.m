function splines = movie2spline(mymovie, type, opts)

  if (nargin < 3)
    fields = {'carth'};
  else
    fields = opts.analysed_fields;
  end

  nframes = length(mymovie.(type).eggshell);
  nfields = length(fields);

  if (isfield(mymovie.(type), 'splines'))
    orig_splines = mymovie.(type).splines;
  else 
    orig_splines = [];
  end
  
  [junk, orig_i, orig_f] = size(orig_splines);
  splines = get_struct('spline',[2 nframes nfields]);

  for i=1:nframes
    for f=1:nfields
      if ((mymovie.(type).update(1,i) | opts.recompute) & isfield(mymovie.(type).eggshell(i), fields{f}))
        splines(1,i,f) = create_spline(mymovie.(type).eggshell(i).(fields{f}));
      elseif (i <= orig_i & f <= orig_f)
        splines(1,i,f) = orig_splines(1,i,f);
      end
      if ((mymovie.(type).update(2,i) | opts.recompute) & isfield(mymovie.(type).cortex(i), fields{f}))
        splines(2,i,f) = create_spline(mymovie.(type).cortex(i).(fields{f}));
      elseif (i <= orig_i & f <= orig_f)
        splines(2,i,f) = orig_splines(2,i,f);
      end
    end
  end

  return;
end
