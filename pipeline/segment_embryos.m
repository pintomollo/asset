function [mymovie, opts] = segment_embryos(myrecording, opts)

  % A nice status-bar if possible
  if (opts.verbosity > 1)
    hwait = waitbar(0,'','Name','CAST');
  end

  % Get the number of channels to parse
  nchannels = length(myrecording.channels);

  % Loop over them
  for indx = 1:nchannels

    % Get the current type of segmentation to apply
    segment_type = myrecording.channels(indx).type;

    % Get the number of frames
    nframes = size_data(myrecording.channels(indx));

    % Prepare the data
    all_embryos = cell(nframes, 1);
    all_estimates = cell(nframes, 1);

    % Update the waitbar
    if (opts.verbosity > 1)
      waitbar(0, hwait, ['Segmenting channel #' num2str(indx) ': ' myrecording.channels(indx).type]);
    end

    % Iterate over the whole recording
    for nimg = 1:nframes

      % We may need data about the noise
      noise = [];

      % Get the current image
      img = double(load_data(myrecording.channels(indx), nimg));

      % Get the noise parameters
      [ellipses, estimates] = perform_step('embryo', segment_type, img, opts);

      % If we have some detections, store them in the final structure
      if (~isempty(ellipses))
        all_embryos{nimg} = ellipses;
        all_estimates{nimg} = estimates;
      end

      % Update the progress bar
      if (opts.verbosity > 1)
        waitbar(nimg/nframes,hwait);
      end
    end

    keyboard

    % Store all detection in the segmentation structure
    new_embryos = filter_embryos(all_embryos, all_estimates, ...
                  opts.split_parameters.angle_thresh, ...
                  opts.split_parameters.max_ratio, ...
                  opts.split_parameters.max_distance / opts.pixel_size, ...
                  opts.split_parameters.max_overlap, ...
                  opts.split_parameters.max_score, opts.split_parameters.max_area_diff);
  end

  % Close the status bar
  if (opts.verbosity > 1)
    close(hwait);
  end


  return;
end
