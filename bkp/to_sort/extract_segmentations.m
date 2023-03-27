function paths = extract_segmentations(mymovie, type, opts)

  if (nargin < 3)
    fields = {'carth'};
  else
    fields = opts.analyzed_fields;
  end

  nframes = length(mymovie.(type).eggshell);
  nfields = length(fields);

  paths = cell([2 nframes nfields]);

  for i=1:nframes
    for f=1:nfields
      if (isfield(mymovie.(type).eggshell(i), fields{f}))
        paths(1,i,f) = {mymovie.(type).eggshell(i).(fields{f})};
      end
      if (isfield(mymovie.(type).cortex(i), fields{f}))
        paths(2,i,f) = {mymovie.(type).cortex(i).(fields{f})};
      end
    end
  end

  return;
end
