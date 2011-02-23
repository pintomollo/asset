function [trackings, expr_name, updated] = import_trackings(trackings, opts)

  updated = false;
  files = opts.trackings;
  expr_name = '';

  %track = struct('mean',[],'name','','errors',[],'expr','','child',[],'files',repmat(struct('fname','','splines',[],'shapes',[],'groups',{}),0,1));
  if (strncmp(opts.segmentation_type, 'all', 3))
    types = {'dic', 'markers'};
  else
    types = {opts.segmentation_type};
  end

  if (isempty(trackings))
    trackings = get_struct('trackings', 1);
  elseif (isstruct(trackings))
    
    done = true;
    for t=1:length(types)
      if (~isfield(trackings, types{t}) & ~isempty(trackings.(types{t})))
        done = false;
      end
    end

    if (done & ~opts.recompute)
      return;
    end
  end

  max_frames = 0;

  for t=1:length(types)
    if (~isempty(files))
      if (ischar(files))
        files = {{files}};
      end

      if (iscell(files))
        for i=1:length(files)
          for j=1:length(files(i))
            if (j==1)
              [tokens,junk]=regexp(files{i}{j},'(.+[-_])?([^-_\.]+)(\..+)','tokens');
              name = tokens{1}{1};
              suffix = tokens{1}{2};
              ext = tokens{1}{3};

              if (strncmp(ext, '.mat', 4))
                tmp = load(fname{f});
                tmp_trackings = tmp.trackings;
                tmp_trackings = tmp_trackings.(types{t});

                for i=1:length(tmp_trackings.child)
                  found = false;
                  for j=1:length(trackings.(types{t}).child)
                    if (strcmp(trackings.(types{t}).child(j).expr, tmp_trackings.child(i).expr))
                      found = true;
                      break;
                    end
                  end

                  if (~found)
                    if (length(trackings.(types{t}).child) == 0)
                      trackings.(types{t}).child = get_struct('tracking');
                    end

                    trackings.(types{t}).child(end+1) = tmp_trackings.child(i);
                  end
                end
              else
                if (isempty(name))
                  expr = [suffix ext];
                  name = suffix;
                else
                  expr = [name '*' ext];
                end

                trackings.(types{t}).child(end+1).name = name;
                trackings.(types{t}).child(end).expr = expr;
                trackings.(types{t}).child(end).files(1).fname = files{i}{j};
              end
            else
              trackings.(types{t}).child(end).files(end+1).fname = files{i}{j};
            end
          end
        end
      elseif (isstruct(files))
        tracking.(types{t}).child = [trackings.(types{t}).child files.(types{t})];
      end
    end

    child_names = {};
    for i=1:length(trackings.(types{t}).child)
      child_names = [child_names trackings.(types{t}).child(i).name];
    end
    [junk, indxs] = sort(child_names);
    trackings.(types{t}).child = trackings.(types{t}).child(indxs);

    if (opts.verbosity > 0)
      expr_name = get_name(trackings.(types{t}), expr_name);

      disp(['[Select your ' expr_name ' ' types{t} ' tracking files]'])

      trackings.(types{t}) = input_trackings(trackings.(types{t}), [types{t} ' (' expr_name ')']);
      expr_name = get_name(trackings.(types{t}), expr_name);
    end
    [trackings.(types{t}), nframes] = load_trackings(trackings.(types{t}), opts);
    if (nframes > max_frames)
      max_frames = nframes;
    end
  end

  for t=1:length(types)
    trackings.(types{t}) = check_trackings(trackings.(types{t}), max_frames);
  end

  updated = true;

  return;
end

function trackings = check_trackings(trackings, max_frames)

  for i=1:length(trackings.files)
    [ngroups,nframes] = size(trackings.files(i).shapes);
    if (nframes < max_frames)
      trackings.files(i).shapes = [trackings.files(i).shapes repmat(struct('path',[]), ngroups, max_frames - nframes)];
      nframes = max_frames;
    end

    for j=1:ngroups
      for k=1:nframes
        if (isempty(trackings.files(i).shapes(j,k).path))
          disp(['File ' trackings.files(i).fname ', tracking of ' trackings.files(i).groups{j} ' in frame ' num2str(k) ' is missing']);
        end
      end
    end
  end

  for i=1:length(trackings.child)
    trackings.child(i) = check_trackings(trackings.child(i), max_frames);
  end

  return;
end

function name = get_name(mystruct, name)

  for i = 1:length(mystruct.child)
    name = get_name(mystruct.child(i), name);
    if (~isempty(name))
      return;
    end
  end

  if (isfield(mystruct, 'name') & ~isempty(mystruct.name))
    indxs = findstr(mystruct.name, '-');
    if (length(indxs) > 1)
      name = mystruct.name(indxs(1)+1:indxs(2)-1);
    end
  end

  return;
end
