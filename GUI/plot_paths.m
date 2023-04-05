function hgroup = plot_paths(h, paths, color)
% PLOT_PATHS draws spot paths as lines.
%
%   HGROUP = PLOT_PATHS(PATHS) draws all PATHS in the current axes, returning a
%   handler to the hggroup HGROUP containing the lines. PATHS should be a cell array
%   each containing one trajectory per spot, ordered as [X_1, Y_1; X_2, Y_2; ...].
%
%   HGROUP = PLOT_PATHS(..., COLORS) defines the color to be used for the lines. A
%   vector/matrix of colors can be provided to differentiate the paths. Default value
%   is 'y'.
%
%   HGROUP = PLOT_PATHS(HGROUP, ...) draws the lines in the provided HGROUP,
%   replacing existing lines (usually faster than creating a new group).
%
%   HGROUP = PLOT_PATHS(HAXES, ...) draws the hggroup in the axes defined by HAXES.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 07.07.2014

  % Input checking and default values
  if (nargin == 1)
    paths = h;
    h = gca;
    color = 'y';
  elseif (nargin == 2)
    if (iscell(h))
      color = paths;
      paths = h;
      h = gca;
    else
      color = 'y';
    end
  end

  % Default value
  if isempty(color)
    color = 'k';
  end

  % For simplicity, we always work with cell arrays
  if (~iscell(paths))
    paths = {paths};
  end

  % Check the number of paths and adapt the number of colors accordingly
  npaths = length(paths);
  if (ischar(color))
    color = [color, color(ones(1, npaths-length(color)))];
    color = color(1:npaths);
  else
    color = [color; color(ones(1, npaths-size(color,1)),:)];
    color = color(1:npaths,:);
  end

  % Check if we got an hggroup directly or if we need to create one
  if (strncmp(get(h, 'Type'), 'hggroup',7))
    hgroup = h;
  else
    hgroup = hggroup('Parent', h);
  end

  % Get the handler to our parent axes
  haxes = get(hgroup, 'Parent');

  % Change the current "hold" status, so we can draw several paths together
  status = get(haxes, 'NextPlot');
  set(haxes,'NextPlot', 'add');

  % Get the handlers and the number of existing paths in the hggroup
  hlines = get(hgroup, 'Children');
  nlines = length(hlines);

  % Need to remember how many lines we have drawn in total
  count = 0;

  % Loog over all paths
  for i = 1:npaths

    % Get the current path
    curr_paths = paths{i};

    % And its color
    if (ischar(color))
      curr_color = color(i);
    else
      curr_color = color(i,:);
    end

    % Draw only visible paths
    if (size(curr_paths, 1) > 1)
      count = count + 1;

      % If we ran out of lines, create a new one
      if (count > nlines)
        line('XData', curr_paths(:,2), 'YData', curr_paths(:,3), 'Parent', hgroup, 'Color', curr_color, 'Marker', 'o');

      % Otherwise, utilise an existing one
      else
        set(hlines(count), 'XData', curr_paths(:,2), 'YData', curr_paths(:,3), 'Color', curr_color);
      end
    end
  end

  % Set back the status
  set(haxes,'NextPlot', status);

  % Delete additional previous circles
  delete(hlines(count+1:nlines))

  % Prevent the output if not needed
  if (nargout == 0)
    clearvars
  end

  return;
end
