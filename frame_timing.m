function times = frame_timing(mymovie, replace_timing)

  if (nargin == 1)

    if (isfield(mymovie, 'dic') & isfield(mymovie.dic, 'metadata'))
      field = 'dic';
    elseif (isfield(mymovie, 'eggshell') & isfield(mymovie.eggshell, 'metadata'))
      field = 'eggshell';
    elseif (isfield(mymovie, 'cortex') & isfield(mymovie.cortex, 'metadata'))
      field = 'cortex';
    elseif (isfield(mymovie, 'data') & isfield(mymovie.data, 'metadata'))
      field = 'data';
    elseif (isfield(mymovie, 'metadata'))
      field = 'tmp';
      mymovie = struct('tmp', mymovie);
    else
      times = [];
      error('None of the expected fields for metadata extraction are present');
    end

    meta = char(mymovie.(field)(1).metadata);

    [tokens] = regexp(meta,'<PlaneTiming DeltaT="(.+?)"', 'tokens');
    times = zeros(1, length(tokens));
    for i=1:length(tokens)
      times(1,i) = str2num(tokens{i}{1});
    end
  else
    ntimings = length(replace_timing);

    fields = fieldnames(mymovie);
    for f = 1:length(fields)
      if (isfield(mymovie.(fields{f}), 'metadata'))
        for j = 1:length(mymovie.(fields{f})) 
          meta = char(mymovie.(fields{f})(j).metadata);
          if (ntimings == 1)
            [tokens] = regexp(meta,'<PlaneTiming DeltaT="(.+?)"', 'tokens');
          end
          [starts, ends] = regexp(meta,'<PlaneTiming DeltaT="(.+?)"');
          for i = length(starts):-1:1
            if (ntimings == 1)
              meta = [meta(1:starts(i)-1) '<PlaneTiming DeltaT="' num2str(str2num(tokens{i}{1}) + replace_timing) '"' meta(ends(i)+1:end)];
            elseif (i <= ntimings)
              meta = [meta(1:starts(i)-1) '<PlaneTiming DeltaT="' num2str(replace_timing(i)) '"' meta(ends(i)+1:end)];
            end
          end
          mymovie.(fields{f})(j).metadata = java.lang.String(meta);
        end
      end
    end


    times = mymovie;
  end

  return;
end
