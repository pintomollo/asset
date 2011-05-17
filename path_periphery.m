function outside_paths = path_periphery(paths, varargin)
% PATH_PERIPHERY computes the periphery of a closed 2D path, thus removing invaginations.
%
% 
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 14.04.2011

  [update, center, axes_length, orientation, opts] = parse_inputs(varargin{:});
  was_cell = true;

  if (~iscell(paths))
    paths = {paths};
    was_cell = false;
  end
  if (isempty(update))
    update = true(size(paths));
  end

  [n,m,o] = size(paths);

  for i=1:n
    for j=1:m
      if (~update(i,j) & ~opts.recompute)
        outside_paths(i,j,:) = paths(i,j,:);

        continue;
      end
      for k=1:o
        if (isempty(paths{i,j,k}))
          continue;
        end

        tmp_path = paths{i,j,k};

        nbins = size(tmp_path, 1);

        if (isempty(center))
          [center, axes_length, orientation] = fit_ellipse(tmp_path);
        end

        outside_theta = [0:nbins]' * (2*pi/nbins);
        outside_pts = zeros(nbins,2);

        polar = carth2elliptic(tmp_path,center,axes_length,orientation); 
        polar = resample(polar, 2);
        %polar = interp_elliptic(polar, full_theta);

        if (all(diff(polar(:,1))>=0))
          % No detectable invagination
          
          outside_paths(i,j,k) = paths(i,j,k);

          continue;
        end

        [junk,ibin] = histc(polar(:,1), outside_theta);

        for b=1:nbins
          bin = polar(ibin==b,:);

          indx = find(bin(:,2)==max(bin(:,2)));
          if (length(indx)>1)
            dist = abs(bin(:,1) - (outside_theta(b) + pi/nbins));
            indx = find(bin(:,2)==max(bin(:,2)) & dist==min(dist),1);
          end
          if (isempty(bin) | isempty(indx))
            outside_pts(b,:) = NaN;
          else
            outside_pts(b,:) = bin(indx,:);
          end
        end

        outside_pts = outside_pts(~any(isnan(outside_pts),2),:);
        outside_pts = elliptic2carth(outside_pts, center, axes_length, orientation);
        outside_paths(i,j,k) = {outside_pts};
      end
    end
  end

  if (~was_cell & numel(paths) == 1)
    outside_paths = outside_paths{1};
  end

  return;
end

function new_pts = resample(pts, factor)

  if (factor < 1)
    new_pts = pts;

    return;
  end

  [npts, ndim] = size(pts);
  new_pts = ones(npts*2, ndim);

  if (abs(pts(1,1)-pts(end,1)) < pi)
    indx = find(pts(1:end-1,1) > pts(2:end,1), 1);
    if (~isempty(indx))
      pts = pts([indx+1:end 1:indx],:);
    end
  end 
  first = pts(1,:);
  first(1,1) = first(1,1) + 2*pi;

  mean_pts = mean(cat(3,pts,[pts(2:end,:); first]), 3);
  new_pts(1:2:end,:) = pts;
  new_pts(2:2:end,:) = mean_pts;

  new_pts = resample(new_pts, factor-1);

  return;
end

function [update, center, axes_length, orientation, opts] = parse_inputs(varargin)

  update = [];
  center = [];
  axes_length = [];
  orientation = [];
  opts = get_struct('ASSET',1);

  if (nargin > 0)
    for i=1:length(varargin)
      type = get_type(varargin{i});
      switch type
        case 'bool'
          update = varargin{i};
        case 'struct'
          opts = varargin{i};
        case 'num'
          if (isempty(center))
            center = varargin{i};
          elseif (isempty(axes_length))
            axes_length = varargin{i};
          elseif (isempty(orientation))
            orientation = varargin{i};
          end
      end
    end
  end

  return;
end
