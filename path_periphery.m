function outside_paths = path_periphery(paths, varargin)
% PATH_PERIPHERY computes the periphery of a closed 2D path, thus removing invaginations.
%
%   OUTSIDE = PATH_PERIPHERY(PATH, OPTS) returns the "outside" part of PATH. This algorithm
%   projects PATH in polar coordinates (see carth2elliptic.m) and returns the outer most
%   point for each angular position, thus getting rid of the invaginations. PATH can
%   either be Nx2 matrix of coordinates or a cell matrix containing multiple paths.
%   Note that this algorithm works only with matrix that have up to 3 dimensions.
%
%   [...] = PATH_PERIPHERY(..., CENTER, AXES_LENGTH, ORIENTATION) projects the paths
%   onto the elliptical coordinate system defined by CENTER, AXES_LENGTH and 
%   ORIENTATION. Note that these three arguments need to be in this particular order
%
%   [...] = PATH_PERIPHERY(..., UPDATE) applies the same algorithm but only for the
%   paths flagged with a TRUE value in the UPDATE matrix. Note that UPDATE must have
%   the same size as the cell matrix PATH in its 2 first dimensions.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 17.05.2011

  % Parse the variable inputs and set the default values
  [update, center, axes_length, orientation, opts] = parse_inputs(varargin{:});
  was_cell = true;

  % We'll work only with cell matrix, so convert it in case
  if (~iscell(paths))
    paths = {paths};
    
    % And remember it so we can return the same output type as the input
    was_cell = false;
  end

  % If it's not defined, we will update all of them
  if (isempty(update))
    update = true(size(paths));
  end

  % Get the 3 first dimensions of the paths
  if (ndims(paths) > 3)
    error('path_pariphery can handle matrices that have up to 3 dimensions only !');
  end
  [n, m, o] = size(paths);

  % Let's loop over all the paths
  for i = 1:n
    for j = 1:m

      % If this path not flagged as to be updated and we do not recompute everything,
      % we are already done !
      if (~any(update(i,j,:)) & ~opts.recompute)
        outside_paths(i,j,:) = paths(i,j,:);

        continue;
      end

      % Loop over the last part
      for k = 1:o
        % If there is nothing in there, do nothing
        if (isempty(paths{i,j,k}))
          continue;
        end

        % Extract the current path
        tmp_path = paths{i,j,k};

        % Get its size
        nbins = size(tmp_path, 1);

        % If we do not have the ellipse provided, estimate it
        if (isempty(center))
          [center, axes_length, orientation] = fit_ellipse(tmp_path);
        end

        % Create the elliptical positions at which we will resample the path
        outside_theta = [0:nbins]' * (2*pi/nbins);
        % And pre-assign the outside points
        outside_pts = zeros(nbins,2);

        % Project the path and over-sample it so that the "bins" contain more than
        % one data point
        polar = carth2elliptic(tmp_path,center,axes_length,orientation); 
        polar = interp_elliptic(polar, 2);

        % No detectable invagination, then continue
        if (all(diff(polar(:,1))>=0))
          outside_paths(i,j,k) = paths(i,j,k);

          continue;
        end

        % We bin the angular position with respect to the defined positions (theta)
        [~, ibin] = histc(polar(:,1), outside_theta);

        % Loop over each bin to get the outter most point
        for b=1:nbins
          % Extract the current points
          bin = polar(ibin == b, :);

          % Get the outter most point in this bin
          indx = find(bin(:, 2) == max(bin(:, 2)));

          % If we have more than one extrema, keep the most central one
          if (length(indx) > 1)
            dist = abs(bin(:, 1) - (outside_theta(b) + (pi / nbins)));
            indx = find((bin(:, 2) == max(bin(:, 2))) & (dist == min(dist)), 1);
          end

          % If we do not have a data point, put a NaN, otherwise store it
          if (isempty(bin) | isempty(indx))
            outside_pts(b,:) = NaN;
          else
            outside_pts(b,:) = bin(indx,:);
          end
        end

        % Keep only the points for which we do have data
        outside_pts = outside_pts(~any(isnan(outside_pts),2),:);

        % Project the new path back onto the cartesian coordinate system
        outside_pts = elliptic2carth(outside_pts, center, axes_length, orientation);

        % Store the result
        outside_paths(i,j,k) = {outside_pts};
      end
    end
  end

  % If we received a matrix, return a matrix
  if (~was_cell & numel(paths) == 1)
    outside_paths = outside_paths{1};
  end

  return;
end

% Sort out the mess due to the variable number of inputs
function [update, center, axes_length, orientation, opts] = parse_inputs(varargin)

  % Initialize the variables with their default values
  update = [];
  center = [];
  axes_length = [];
  orientation = [];
  opts = get_struct('ASSET',1);

  % Loop over the variable inputs
  if (nargin > 0)
    for i = 1:length(varargin)

      % Get the type of the current input and assign it accordingly
      type = get_type(varargin{i});
      switch type

        % A boolean is for the update matrix
        case 'bool'
          update = varargin{i};

        % A structure is the option
        case 'struct'
          opts = varargin{i};

        % For numerical values, we assume they are provided in the correct order:
        % center, axes_length and orientation
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
