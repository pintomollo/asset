function [new_angles, path] = interp_elliptic(orig_angles, orig_values, new_angles)
% INTERP_ELLIPTIC linear interpolation of the values in elliptical coordinate system,
% taking into account the fact that the angular new_anglesitions are periodic.
%
%   [VALUES] = INTERP_ELLIPTIC(ORIGINAL, NEW_POSITION) interpolates ORIGINAL at
%   NEW_POSITION and returns the result in VALUES. Both VALUES and ORIGINAL are NxM
%   matrices in which the first column is the angular new_anglesition.
%
%   [...] = INTERP_ELLIPTIC(ORIGINAL, NRESAMPLE) resamples the original values with
%   a seemingly 2^NRESAMPLE higher resolution.
%
%   [ANGULAR_POSITION, VALUES] = INTERP_ELLIPTIC(...) returns the angular new_anglesition
%   and the interpolated values separatly.
%
%   [...] = INTERP_ELLIPTIC(ORIG_POSITION, ORIG_VALUES, NEW_POSITION) provides the
%   original angular new_anglesitions and values separatly.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 18.05.2011

  % Handle the different types of new_anglessible inputs
  % For 2, we have both the new_anglesition and the values together
  if (nargin == 2)
    new_angles = orig_values;

    orig_values = orig_angles(:, 2:end);
    orig_angles = orig_angles(:, 1);

  % One argument is not valid, so return directly
  elseif (nargin == 1)
    new_angles = orig_angles;

    return;
  end

  % Vectorize the input angles
  new_angles = new_angles(:);

  % Initialization
  is_resampling = false;

  % If we have only one value for the angles, this is a resampling problem
  if (numel(new_angles) == 1 & mod(new_angles, 1) == 0)
    is_resampling = true;

    % Re-aling the positions as the interpolation works only with monotonically
    % increasing positions
    resample_indx = find(orig_angles(1:end-1) > orig_angles(2:end), 1);
    if (~isempty(resample_indx))
      orig_angles = orig_angles([resample_indx+1:end 1:resample_indx]);
      orig_values = orig_values([resample_indx+1:end 1:resample_indx], :);
    end

    % Prepare the temporary variables that we need to re-sample
    resample = new_angles;
    new_angles = orig_angles;

    % At each loop we'll increase the sampling by 2
    for i = 1:resample

      % Create the new vector
      tmp_new_angles = NaN(size(new_angles,1)*2, 1);

      % Compute the very last element, to get periodicity
      first = new_angles(1);
      first = first + 2*pi;

      % Simply average the following positions to get the newer intermediater one
      mean_new_angles = mean([new_angles,[new_angles(2:end,:); first]], 2);

      % Store them
      tmp_new_angles(1:2:end-1,:) = new_angles;
      tmp_new_angles(2:2:end,:) = mean_new_angles;

      % Use this new position vector as the interpolation "target"
      new_angles = tmp_new_angles;
    end
  end

  % Measure the average angle difference between neighboring points
  dangle = mean(diff(new_angles(:,1))); 

  % In case we have a negative value, there is a major problem
  if (any(dangle) < 0)
    error('Interpolation can only be performed on monotonically increasing positions')
  end

  % Deduce the minimal and maximal positions we can interpolate
  min_angle = min(new_angles(:,1)) - dangle;
  max_angle = max(new_angles(:,1)) + dangle;

  % Reset the angular values such that they are in the interpolated range
  orig_angles(orig_angles < min_angle) = orig_angles(orig_angles < min_angle) + 2*pi;
  orig_angles(orig_angles > max_angle) = orig_angles(orig_angles > max_angle) - 2*pi;

  % Re-align the data such that they ahave increasing angular values
  indx = find(orig_angles(1:end-1) > orig_angles(2:end), 1);
  if (~isempty(indx))
    orig_angles = orig_angles([indx+1:end 1:indx]);
    orig_values = orig_values([indx+1:end 1:indx], :);
  end

  % Get the size of the problem
  npts = length(orig_angles);

  % Get the number of points we need to repeat on both sides to avoid border effects
  start = find(orig_angles - 2*pi < new_angles(1,1), 1, 'last');
  ends = find(orig_angles + 2*pi > new_angles(end,1), 1, 'first');

  % Repeat the data on both sides to avoid border effects when interpolating
  orig_angles = orig_angles([start:end 1:end 1:ends]);
  orig_values = orig_values([start:end 1:end 1:ends], :);
  orig_angles(1:npts-start+1) = orig_angles(1:npts-start+1) - 2*pi;
  orig_angles(end-ends+1:end) = orig_angles(end-ends+1:end) + 2*pi;

  % Perform the actual linear interpolation
  path = interp1q(orig_angles, orig_values, new_angles);

  % If we did some re-sampling and had to re-align the vector, switch back to its
  % original alignment
  if (is_resampling && ~isempty(resample_indx))
    resample_indx = resample_indx * (2^new_angle);

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
