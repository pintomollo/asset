function [trackings, max_frames] = load_trackings(trackings, opts)
% LOAD_TRACKINGS loads the content of the .shapes files into the provided structure.
%
%   TRACKINGS = LOAD_TRACKINGS(TRACKINGS, OPTS) returns the same TRACKING structure
%   to which the content of the .shapes files will be added.
%
%   [TRACKINGS, NFRAMES] = LOAD_TRACKINGS(...) also returns the maximum number of
%   frames among all the .shapes files.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 23.03.2011

  % Initialize the output variable
  max_frames = 0;

  % If trackings has not the correct structure, try to find a valid sub-structure
  if (~isfield(trackings, 'child'))
    % List the fields
    fields = fieldnames(trackings);

    % Loop over the fields
    for i=1:length(fields)
      field = fields{i};

      % If it's a structure, load it recursively
      if (isstruct(trackings.(field)))
        trackings.(field) = load_trackings(trackings.(field), opts);
      end
    end

    return;
  end

  % If a 'shapes' field exists, we are the the correct level of the structure
  if (isfield(trackings.child, 'shapes'))

    % Loop over all the children of the current group
    for i=1:length(trackings.child)

      % If we need to reload (or if there is no data yet)
      if (isempty(trackings.child(i).shapes) || opts.recompute)
        % Make sure the path is relative
        trackings.child(i).fname = relativepath(trackings.child(i).fname);

        % Load the .shape file
        [trackings.child(i).shapes, trackings.child(i).groups] = load_shapes(trackings.child(i).fname);
      end

      % If we need to recompute the periphery, do it
      if (opts.recompute & opts.follow_periphery)
        trackings.child(i).shapes = path_periphery(trackings.child(i).shapes);
      end

      % Compute the number of frames and update the max_frames if needed
      nframes = size(trackings.child(i).shapes,2);
      if (nframes > max_frames)
        max_frames = nframes;
      end
    end

  % Otherwise, we are not yet at the correct level
  else
    % We have a look in every child recursively
    for i=1:length(trackings.child)
      [trackings.child(i), nframes] = load_trackings(trackings.child(i), opts);

      % Update the number of frames if needed
      if (nframes > max_frames)
        max_frames = nframes;
      end
    end
  end

  % Update the expression
  if (isfield(trackings, 'expr'))
    trackings.expr = relativepath(trackings.expr);
  end

  return;
end
