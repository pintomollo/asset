function colors = colorize_graph(coords, colors)
% COLORIZE_GRAPH tries to assign colors as different as possible to neighboring vertices
% of a graph, using a very simple heuristic.
%
%   GCOLORS = COLORIZE_GRAPH(COORDS, COLORS) assign one COLORS to each COORDS, intenting
%   to avoid similar COLORS on neighboring COORDS. GCOLORS has as many rows as COORDS.
%
%   GCOLORS = COLORIZE_GRAPH(COORDS) utilizes 'redbluemap' to obtain the COLORS to
%   choose from.
%
%   GCOLORS = COLORIZE_GRAPH(PATHS, ...) computes the average position of the various
%   PATHS to compute their distance to neighbors.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 07.07.2014

  % Double check that we have some data
  if (nargin == 0 || isempty(coords))
    colors = [];

    return;
  end

  % Check if we need to get colors
  if (nargin == 1)
    if (iscell(coords))
      colors = redbluemap(length(coords));
    else
      colors = redbluemap(size(coords, 1));
    end
  end

  % Compute the average position for each path
  if (iscell(coords))
    paths = coords;

    empty_cells = cellfun('isempty', paths);

    coords = cellfun(@(x)(mymean(x(:,2:3), 1)), paths(~empty_cells), 'UniformOutput', false);
    coords = cat(1, coords{:});
  end

  % We do not have anything to work on
  if (isempty(coords))
    return;
  end

  % Verify that we do not have doublets (because of the integer pixel values)
  [tmp_coord, junk, indx] = unique(coords, 'rows');
  if (numel(tmp_coord) ~= numel(coords))

    % In this case, we loop over the duplicates and addsome random noise
    for i=1:length(junk)
      goods = (indx==i);
      if (sum(goods)>1)
        coords(goods, :) = coords(goods, :) + rand(sum(goods), 2);
      end
    end
  end

  % Split the coordinates
  xcoord = coords(:,1);
  ycoord = coords(:,2);

  % Compute their mutual distance
  dist = (bsxfun(@minus, xcoord, xcoord.').^2 + bsxfun(@minus, ycoord, ycoord.').^2);

  icolors = [];
  % Try to simplify the all-to-all distance, using Delaunay triangulation
  if (length(xcoord) > 2)
    try
      trig = delaunay(xcoord, ycoord);

      % Get the respective matrix coordinates to retrieve the corresponding distances
      icoord = trig;
      jcoord = trig(:,[2:end 1]);

      icoord = icoord(:);
      jcoord = jcoord(:);

      % Sort the distances
      subs = sub2ind(size(dist), icoord, jcoord);
      edge_dists = dist(subs);
      [junk, indx] = sort(edge_dists);

      % Reorder the indexes
      icoord = icoord(indx);
      jcoord = jcoord(indx);

      % Get the number of colors and the two indexes used to choose a new color
      ncolors = size(colors, 1);
      sindx = 1;
      eindx = floor(ncolors/2)+1;

      % Initialize the index array for colors
      icolors = NaN(size(xcoord));

      % Loop over all edges
      for i=1:length(icoord)

        % Anything left to do ?
        todo = isnan(icolors);
        invert = false;

        % Assign the current vertex ?
        if (todo(icoord(i)))

          % Just check if it is closer to the other index
          subdists = dist(icoord(i), todo);
          [v, indx] = min(subdists);

          % If there has been chosen colors, find the closest one
          if (~isempty(indx))
            tmpc = icolors(todo);
            tmpc = tmpc(indx);

            % Maybe we better choose the other index
            invert = (abs(tmpc - sindx) < abs(tmpc - eindx));
          end

          % Assign a color from one of the two index, and update them
          if (invert)
            icolors(icoord(i)) = eindx;
            eindx = mod(eindx, ncolors)+1;
          else
            icolors(icoord(i)) = sindx;
            sindx = mod(sindx, ncolors)+1;
          end
        end

        % If we need to find a color, choose it from one of the two indexes
        if (todo(jcoord(i)))
          if (invert)
            icolors(jcoord(i)) = sindx;
            sindx = mod(sindx, ncolors)+1;
          else
            icolors(jcoord(i)) = eindx;
            eindx = mod(eindx, ncolors)+1;
          end
        end

        % All were chosen
        if (~any(todo))
          break;
        end
      end
    catch
      % Nothing
    end
  end

  % Something went wrong...
  if (isempty(icolors))
    icolors = randperm(sum(~empty_cells)).';
  end

  % Reorder colors and return it
  colors(~empty_cells, :) = colors(icolors, :);

  return;
end
