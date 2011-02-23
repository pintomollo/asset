function errors = path_error(ref, paths, varargin)

  if (numel(paths) == 0)
    errors = NaN;

    return;
  end

  [update, orig_errors, centers, axes_length, orientations, nbins, opts] = parse_inputs(varargin{:});

  [ntypes, nframes, npaths] = size(paths);
  if (isempty(update))
    update = true(size(paths));
  end
  if (isempty(orig_errors))
    orig_errors = NaN(size(paths));
  end

  theta = (2*pi) / nbins;
  bins = [0 0.5:nbins nbins] * theta;
  errors = NaN(ntypes, nframes, npaths, nbins);
  
  for t=1:ntypes
    if (~any(update(t,:)) & ~opts.recompute)
      errors(t,:,:,:) = orig_errors(t,:,:,:);

      continue;
    end

    for f=1:nframes
      if (~update(t,f) & ~opts.recompute)
        errors(t,f,:,:) = orig_errors(t,f,:,:);
  
        continue;
      end

      tmp_ref = ref{t,f};
  
      if (~isempty(tmp_ref))
      %  errors(t,f,:,:) = NaN;
      %else
        if (any(isnan(centers(:,f))))
          [centers(:,f), axes_length(:,f), orientations(1,f)] = fit_ellipse(tmp_ref);
        end
        ell_ref = carth2elliptic(tmp_ref, centers(:,f), axes_length(:, f), orientations(1,f));

        for i=1:npaths
          if (isempty(paths{t,f,i}))
            errors(t,f,i,:) = NaN;
          else
            ell_path = carth2elliptic(paths{t,f,i}, centers(:,f), axes_length(:, f), orientations(1,f));
            areas = compute_error(ell_ref, ell_path);

            errors(t,f,i,:) = bin_error(areas, bins);
          end
        end
      end
    end
  end
  
  errors = errors / (2*pi);

  return;
end

function areas = compute_error(ref, path)

  path = interp_elliptic(path, ref(:,1));
  dt = diff(ref(:,1));

  dr = path(:,2) - ref(:,2);
  ds = abs(diff(sign(dr)));
  dr = abs(dr);

  areas = (dr(1:end-1) + dr(2:end)) .* dt ./ 2;

  if (any(ds ~= 0))
    targets = (ds ~= 0);
    denum = dt(targets) ./ (2* (dr([targets; false]) + dr([false; targets])));
    indxs = ~isinf(denum);

    denum = denum(indxs);
    targets(~indxs) = false;
    
    areas(targets) = (dr([targets; false]).^2 + dr([false; targets]).^2) .* denum;
  end
  areas = [ref(1:end-1,1)+dt abs(areas)];

  return;
end

function errors = bin_error(areas, bins)

  nbins = length(bins) - 2; 
  [junk, map] = histc(areas(:,1), bins);
  errors = zeros(1,nbins);

  for k=1:nbins
    errors(k) = sum(areas(map == k, 2));
  end
  errors(1) = errors(1) + sum(areas(map(1:end-1) == (nbins + 1), 2));

  return;
end

function [update, orig_errors, centers, axes_length, orientations, nbins, opts] = parse_inputs(varargin)

  update = [];
  orig_errors = [];
  centers = [];
  axes_length = [];
  orientations = [];
  opts = get_struct('ASSET', 1);
  nbins = [];

  if (nargin > 0)
    for i=1:length(varargin)
      type = get_type(varargin{i});
      switch type
        case 'bool'
          update = varargin{i};
        case 'num'
          if (numel > 1)
            orig_errors = varargin{i};
          else
            nbins = varargin{i};
          end
        case 'struct'
          if (isfield(varargin{i}, 'centers'))
            centers = varargin{i}.centers;
            axes_length = varargin{i}.axes_length;
            orientations = varargin{i}.orientations;
          else
            opts = varargin{i};
          end
      end
    end
  end

  if (isempty(nbins))
    nbins = opts.nbins;
  end

  return;
end
