function [new_angles, path] = interp_elliptic(varargin)
% INTERP_ELLIPTIC linear interpolation of the values in elliptical coordinate system,
% taking into account the fact that the angular positions are periodic.
%
%   [VALUES] = INTERP_ELLIPTIC(ORIGINAL, NEW_POSITION) interpolates ORIGINAL at
%   NEW_POSITION and returns the result in VALUES. Both VALUES and ORIGINAL are NxM
%   matrices in which the first column is the angular position.
%
%   [...] = INTERP_ELLIPTIC(ORIGINAL, NRESAMPLE, IS_RESAMPLING) resamples the original
%   values with a seemingly 2^NRESAMPLE higher resolution when IS_RESAMPLING is true.
%
%   [ANGULAR_POSITION, VALUES] = INTERP_ELLIPTIC(...) returns the angular position
%   and the interpolated values separatly.
%
%   [...] = INTERP_ELLIPTIC(ORIG_POSITION, ORIG_VALUES, ...) provides the original 
%   angular positions and values separatly.
%
%   [...] = INTERP_ELLIPTIC(..., ANGLE_RANGE) interpolates in between ANGLE_RANGE
%   instead of the default [0 2*pi] range.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 18.05.2011

  [orig_angles, orig_values, new_angles, angle_range, is_resampling] = parse_input(varargin{:});

  % Measure the "size" of the range
  drange = diff(angle_range);

  % Vectorize the input angles
  new_angles = new_angles(:);

  % If we have only one value for the angles, this is a resampling problem
  if (is_resampling)

    % Re-aling the positions as the interpolation works only with monotonically
    % increasing positions
    resample_indx = find(orig_angles(1:end-1) - orig_angles(2:end) > 0.5*drange, 1);
    if (~isempty(resample_indx))
      orig_angles = orig_angles([resample_indx+1:end 1:resample_indx]);
      orig_values = orig_values([resample_indx+1:end 1:resample_indx], :);
    end

    % Prepare the temporary variables that we need to re-sample
    resample = new_angles;
    new_angles = orig_angles;
    new_values = orig_values;

    % At each loop we'll increase the sampling by 2
    for i = 1:resample

      % Create the new vector
      tmp_new_angles = NaN(size(new_angles,1)*2, 1);
      tmp_new_values = NaN(size(new_values,1)*2, size(new_values, 2));

      % Compute the very last element, to get periodicity
      first = new_angles(1);
      first = first + drange;
      first_val = new_values(1,:);

      % Simply average the following values to get the newer intermediate ones
      mean_new_angles = mean([new_angles,[new_angles(2:end,:); first]], 2);
      mean_new_values = mean(cat(3,new_values,[new_values(2:end,:); first_val]), 3);

      % Store them
      tmp_new_angles(1:2:end-1,:) = new_angles;
      tmp_new_angles(2:2:end,:) = mean_new_angles;
      tmp_new_values(1:2:end-1,:) = new_values;
      tmp_new_values(2:2:end,:) = mean_new_values;

      % Use this new position vector as the interpolation "target"
      new_angles = tmp_new_angles;
      new_values = tmp_new_values;
    end
    path = new_values;

    % If we have only one output, combine both
    if (nargout == 1)
      new_angles = [new_angles path];
      path = [];
    end

    return;
  end

  % Measure the average angle difference between neighboring points
  dangle = mean(diff(new_angles(:,1))); 

  % In case we have a negative value, there is a major problem
  if (any(dangle) < 0)
    error('Interpolation can only be performed on monotonically increasing positions')
  end

  % Deduce the minimal and maximal positions we can interpolate
  min_angle = angle_range(1) - dangle;
  max_angle = angle_range(2) + dangle;

  % Reset the angular values such that they are in the interpolated range
  orig_angles(orig_angles < min_angle) = orig_angles(orig_angles < min_angle) + drange;
  orig_angles(orig_angles > max_angle) = orig_angles(orig_angles > max_angle) - drange;

  % Remove the small backward values
  smalls = [true false];
  while (any(smalls) & any(~smalls))
    dorig = diff(orig_angles([1:end 1]));
    smalls = (dorig <= 0 & dorig > -drange/4);
    orig_angles = orig_angles(~smalls);
    orig_values = orig_values(~smalls, :);
  end

  if (isempty(orig_angles))
    % If we have only one output, combine both
    if (nargout == 1)
      new_angles = [new_angles NaN(size(new_angles))];
    end

    return;
  end

  % Re-align the data such that they have increasing angular values
  indx = find(orig_angles(1:end-1) > orig_angles(2:end));

  if (length(indx) > 1)
    warning('Non-monotonic angular values, interpolation can produce incoherent values');
  elseif (~isempty(indx))
    orig_angles = orig_angles([indx+1:end 1:indx]);
    orig_values = orig_values([indx+1:end 1:indx], :);
  end

  % Get the size of the problem
  npts = length(orig_angles);

  % Get the number of points we need to repeat on both sides to avoid border effects
  start = find(orig_angles - drange < new_angles(1,1), 1, 'last');
  ends = find(orig_angles + drange > new_angles(end,1), 1, 'first');

  % Repeat the data on both sides to avoid border effects when interpolating
  orig_angles = orig_angles([start:end 1:end 1:ends]);
  orig_values = orig_values([start:end 1:end 1:ends], :);
  orig_angles(1:npts-start+1) = orig_angles(1:npts-start+1) - drange;
  orig_angles(end-ends+1:end) = orig_angles(end-ends+1:end) + drange;

  % Perform the actual linear interpolation
  path = interp1q(orig_angles, orig_values, new_angles);

  % If we did some re-sampling and had to re-align the vector, switch back to its
  % original alignment
  if (is_resampling && ~isempty(resample_indx))
    resample_indx = resample_indx * (2^resample);

    path = path([end-resample_indx+1:end 1:end-resample_indx], :);
    new_angles = new_angles([end-resample_indx+1:end 1:end-resample_indx], :);
  end

  % If we have only one output, combine both
  if (nargout == 1)
    new_angles = [new_angles path];
    path = [];
  end

  return;
end

% Sort out the mess due to the variable number of inputs
function [orig_angles, orig_values, new_angles, angle_range, is_resampling] = parse_input(varargin)

  % Initialize the output variable
  orig_angles = NaN;
  orig_values = NaN;
  new_angles = NaN;
  angle_range = NaN;
  is_resampling = false;

  % Loop over all the input parameters
  for i=1:length(varargin)

    % We assign them in a specific order, checking whether they've been assigned,
    % first the original values, then the new positions
    if (isnan(orig_values))
      orig_values = varargin{i};
    elseif (isnan(new_angles))
      new_angles = varargin{i};

    % Here it becomes more complicated as this can mean several things !
    else
      % A boolean has to be the resampling flag
      if (islogical(varargin{i}))
        is_resampling = varargin{i};

      % A short vector is the angle range
      elseif (numel(varargin{i}) == 2 & size(orig_values, 1) ~= size(new_angles, 1))
        angle_range = varargin{i};

      % Otherwise it has to be the new positions while the two previous vectors
      % were the values and the angles separate
      else
        orig_angles = orig_values;
        orig_values = new_angles;
        new_angles = varargin{i};
      end
    end
  end

  % If we still do not have the angles, they are in the values
  if (isnan(orig_angles))
    orig_angles = orig_values(:, 1);
    orig_values = orig_values(:, 2:end);
  end

  % The default range is the entire ellipse
  if (isnan(angle_range))
    angle_range = [0 2*pi];
  end

  return;
end
