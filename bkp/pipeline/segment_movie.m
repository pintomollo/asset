function [mymovie, updated] = segment_movie(mymovie, opts)
% SEGMENT_MOVIE performs the segmentation of the provided recording, detecting the
% position of the eggshell and the cellular membrane in each frame independently.
% By default, this function will try to avoid recomputing the segmentation. To force
% it, set the parameter 'recompute' to true (help load_parameters.m).
%
%   MYMOVIE = SEGMENT_MOVIE(MYMOVIE, OPTS) segments MYMOVIE using the parameters
%   specified in OPTS (see get_struct('ASSET')). The results are stored in MYMOVIE in
%   the field corresponding to the segmented channel.
%
%   [MYMOVIE, UPDATED] = SEGMENT_MOVIE(...) returns whether MYMOVIE has been UPDATED
%   by this function.
%
%   [...] = SEGMENT_MOVIE(MYMOVIE) uses the default value for OPTS as provided by
%   get_struct('ASSET').
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 02.09.2011

  % Nothing has changed yet
  updated = false;

  % Set the default value if need be
  if (nargin < 2)
    opts = get_struct('options', 1);
  end

  % Get the size of the problem, using the correct channel
  switch (opts.segmentation_type)

    % The default channel which we should always have
    case {'dic', 'all'}
      [nframes imgsize ] = size_data(mymovie.dic);

    % The marker segmentation is a bit more problematic as there might be no channel
    % for the eggshell
    case 'markers'
      if (isfield(mymovie, 'eggshell') & ~isempty(mymovie.eggshell))
        [nframes imgsize ] = size_data(mymovie.eggshell);
      else
        [nframes imgsize ] = size_data(mymovie.cortex);
      end

    case 'data'
      [nframes, imgsize] = size_data(mymovie.data);

    % Somehow the input parameters do not make sense
    otherwise
      error 'None of the expected field are present in ''mymovie''';
  end

  % Try as hard as possible not to recompute the segmentation !!
  % Basically if there is some data already, we assume they are good and keep them.
  if (~opts.recompute && strncmp(opts.do_ml, 'none', 4)) && ...
     (((strncmp(opts.segmentation_type, 'dic', 3) | strncmp(opts.segmentation_type, 'all', 3)) && ...
       (isfield(mymovie.dic, 'eggshell') && isfield(mymovie.dic, 'cortex')) && ...
       (length(mymovie.dic.eggshell) == nframes && length(mymovie.dic.cortex) == nframes)) ...
     || ...
     (strncmp(opts.segmentation_type, 'data', 3) && ...
       (isfield(mymovie.data, 'eggshell') && isfield(mymovie.data, 'cortex')) && ...
       (length(mymovie.data.eggshell) == nframes && length(mymovie.data.cortex) == nframes)) ...
     || ...
     ((strncmp(opts.segmentation_type, 'markers', 7) | strncmp(opts.segmentation_type, 'all', 3)) && ...
       (isfield(mymovie.markers, 'eggshell') && isfield(mymovie.markers, 'cortex')) && ... 
       (length(mymovie.markers.eggshell) == nframes && length(mymovie.markers.cortex) == nframes)))

    return;
  end

  % Progress bar
  if (opts.verbosity > 0)
    hwait = waitbar(0,'Segmenting Frames','Name','ASSET');
  end

  switch (opts.parse_frames)
    case 'normal'
      frames = 1:nframes; 
    case 'random'
      frames = randperm(nframes);      
  end

  if (opts.max_frames > 1)
    max_frames = round(opts.max_frames);
  else
    max_frames = round(nframes*opts.max_frames);
  end

  if max_frames > nframes 
    max_frames = nframes;
  elseif max_frames < 1
    max_frames = 1;
  end

  frames = frames(1:max_frames);
  dic_opts = opts;
  dic_opts.segmentation_type = 'dic';

  if (opts.recompute)
    mymovie = split_cells(mymovie, opts);
  end

  did_dic = false;
  did_data = false;
  for i = 1:max_frames
    nframe = frames(i);

    if (strncmp(opts.segmentation_type, 'dic', 3) | strncmp(opts.segmentation_type, 'all', 3))
      mymovie = dp_dic(mymovie, nframe, opts);

      updated = updated || any(mymovie.dic.update(:, nframe));
    end

    if (strncmp(opts.segmentation_type, 'data', 3))
      did_data = true;
      mymovie = dp_data(mymovie, nframe, opts);

      updated = updated || any(mymovie.data.update(:, nframe));
    end

    if (strncmp(opts.segmentation_type, 'markers', 7) | strncmp(opts.segmentation_type, 'all', 3))
      if (opts.recompute |~isfield(mymovie, 'eggshell') | isempty(mymovie.eggshell))
        did_dic = true;
        mymovie = dp_dic(mymovie, nframe, opts);

        if (any(mymovie.dic.update(:, nframe)))

          tmp_opts = dic_opts;
          tmp_opts.recompute = true;
          mymovie = duplicate_segmentation(mymovie, 'markers', tmp_opts, nframe);
          tmp_opts = opts;
          tmp_opts.recompute = true;

          mymovie = dp_markers(mymovie, nframe, tmp_opts);
        else
          mymovie.markers.update = mymovie.dic.update;
        end
      else
        mymovie = dp_markers(mymovie, nframe, opts);
      end

      updated = updated || any(mymovie.markers.update(:, nframe));
    end

    if (opts.verbosity > 0)
      waitbar(i/max_frames,hwait);
    end
  end

  if (did_dic)
    mymovie = correct_jitter(mymovie, dic_opts);
    mymovie = align_embryo(mymovie, dic_opts);
  end
  mymovie = correct_jitter(mymovie, opts);
  mymovie = align_embryo(mymovie, opts);

  if (opts.verbosity > 0)
    close(hwait);
  end

  return;
end
