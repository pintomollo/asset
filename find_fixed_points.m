function [fixed_points] = find_fixed_points(Fiso, range, function_params)
% FIND_FIXED_POINTS studies the fixed_points of 2D differential system in the provided 
% range and returns their position.
%
%   [FPS] = FIND_FIXED_PONTS(FISO, RANGE, PARAMS) returns the position of the fixed
%   points found in RANGE as a Nx2 matrix in FPS. To do so, it uses the provided
%   FISO which is a function handler for the two isoclines, along with the provided
%   function parameters PARAMS.
%
%     FISO              Function handler of the isoclines of the differential 
%                       equation: dx/dt = 0 -> y = FISO(x, PARAMS)
%                       They should be organized such as to return the first
%                       isocline, then the second one, in the 3rd dimension
%                       (i.e. cat(3, [isox_x, isox_y], [isoy_x, isoy_y]))
%

  % Initialize the output variables to avoid runtime errors
  fixed_points = [];

  % Some functions do not have parameters so an empty vector will do
  if (nargin < 3)
    function_params = [];
  end

  % We need a valid range to work with
  if (numel(range) < 2 | range(2) <= range(1))
    disp('A valid interval needs to be provided !');

    return;
  end

  % If we have 4 elements in range, we have two asymmetric dimensions
  if (numel(range) == 4)
    % Discretize the interval in the two dimensions
    pts = [linspace(range(1), range(2)).' linspace(range(3), range(4)).'];
  else
    % Use a symmetric discretization
    pts = linspace(range(1), range(2)).';
    pts = pts(:,[1 1]);
  end

  % Now let's estimate the position of the fixed points (FP) !
  % We check where the difference between the two isoclines change (help sign)
  sign_iso = sign(cross(Fiso, pts, function_params));

  % When the sign change, there must be a FP in between (help diff)
  cross_pos = abs(diff(sign_iso)); 

  % Now we look for intervals in between which the sign changes (i.e. the difference
  % is not 0 !)
  pts_pos = find(cross_pos > 0);

  % Count how many FP we found
  npts = length(pts_pos);

  % Prepare the output matrix (nptsx2)
  fixed_points = NaN(npts, 2);

  % Find each FP separatly
  for i=1:npts

    % Index of the interval in which it is located
    index = pts_pos(i);

    % If the sign is 0, we are at the fixed point
    if (sign_iso(index) == 0)
      tmp_pos = pts(index, :);

    % If the next point is exactly at the fixed point, we'll handle this later
    elseif (sign_iso(index+1) == 0)
      continue;
    else
      % Finally, we use fzero (help fzero) to precisely find the position of the FP.
      % To do this, we use the fact that cross will be 0 at the FP !
      tmp_pos = fzero(@(u) cross(Fiso, u, function_params), pts([index index+1]));
    end

    % Now that we know its position in one dimension, we use the isocline
    % to find the other coordinate and store it. It does not matter from which
    % isocline we copy the position as, by definition, both are the same
    tmp_pos = Fiso(tmp_pos(1,[1 1]), function_params);
    fixed_points(i, :) = tmp_pos(1, :, 1);
  end

  % Remove possible duplicates
  fixed_points = fixed_points(any(~isnan(fixed_points), 2), :);

  return;
end

% Help function that computes the difference between the two isoclines
function dist = cross(Fiso, x, function_params)

  % We need a two-colon matrix, thus ensure this !
  if (size(x, 2) == 1)
    x = [x x];
  end

  % Compute the value of the isocline
  pts = Fiso(x, function_params);

  % We need to increase its dimensionalty to be able to compare both isoclines
  x = cat(3, x, x);

  % Identify which dimension has been modified (i.e. y = f(x) or x = f(y) ?)
  [index, junk] = find(squeeze(any(x ~= pts, 1)));

  % If we have already one isocline ready, focus on the second one
  if (numel(index) == 1)
    index = index([1 1]);
  end

  % No point different, we are at the fixed point already
  if (isempty(index))
    % Thus distnce is null
    dist = zeros(size(x, 1), 1);

  % Both isocline work with the same coordinate
  elseif (index(1) == index(2))
    
    % Simply substract one from the other
    dist = pts(:, index(1), 1) - pts(:, index(2), 2);

  % They use different coordinates, we need to feed the other isocline with the result
  % of the first one
  else

    % Find out which isocline has re-computed the first coordinate
    [junk, isoclines] = sort(index);

    % Compute the reverse position (i.e. the position as estimated by the other isocline)
    reverse_pts = Fiso([pts(:, 1, isoclines(1)) pts(:, 2, isoclines(2))], function_params);

    % Measure the distance with the initial coordinate
    dist = x(:, index(2), 1) - reverse_pts(:, index(2), 2);
  end

  return;
end
