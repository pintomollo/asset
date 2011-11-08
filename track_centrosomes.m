function mymovie = track_centrosomes(mymovie, spots, opts)
% TRACK_CENTROSOMES detects and tracks over time the position of the two centrosomes
% in a fluorescence time-lapse recording. This algorithm is designed to find two and
% only two "particles".
%
%   MYMOVIE = TRACK_CENTROSOMES(MYMOVIE, OPTS) detects particles in the "data"
%   channel of MYMOVIE, identifies the centrosomes and tracks them over time. These
%   results are stored in the "centrosomes" field of "data".
%
%   [...] = TRACK_CENTROSOMES(MYMOVIE) uses the default value of OPTS as provided by
%   get_struct('ASSET').
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 01.09.2011

  if (nargin == 2)
    opts = spots;
    spots = {{}};
  end

  % Check whether we actually do have to compute anything 
  if (isfield(mymovie.data, 'centrosomes') & ~isempty(mymovie.data.centrosomes) & ~opts.recompute)
    return;
  end

  % Make sure the data have been copied to the current field
  if (~isfield(mymovie.data, 'centers') | isempty(mymovie.data.centers))
    mymovie = duplicate_segmentation(mymovie, 'data', opts);
    %save(mymovie.experiment, 'mymovie', 'trackings','opts');
  end

  if (isfield(mymovie.data, 'projection') & ~isempty(mymovie.data.projection) & exist(mymovie.data.projection.fname, 'file') == 2)

    [nframes, imgsize] = size_data(mymovie.data.projection);

    if (isempty(spots))
      spots = detect_spots(mymovie.data.projection, opts);
    end
  else
    % The recording size
    [nframes, imgsize] = size_data(mymovie.data);
    % Identify the candidate spots
    if (isempty(spots))
      spots = detect_spots(mymovie.data, opts);
    end
  end

  %%%%%%%%% SHOULD NOT BE HARD-CODED
  volume_thresh = 0.01;
  score_thresh = 0.001;

  % We are looking for only a subset of the detected spots, let's filter them
  for i = 1:nframes
    spot_volume = 2*pi*spots{i}(:,4).*(spots{i}(:,3).^2);
    spot_score = spots{i}(:,6);

    good_spots = spots{i}(spot_volume > volume_thresh, :);
    better_spots = spots{i}(spot_volume > volume_thresh & spot_score > score_thresh, :);

    if (false)
      imshow(imnorm(double(load_data(mymovie.data.projection, i))));
      hold on;
      scatter(spots{i}(:,2), spots{i}(:, 1), 'b');
      scatter(good_spots(:,2), good_spots(:, 1), 'r');
      scatter(better_spots(:,2), better_spots(:, 1), 'y');
      hold off;

      saveas(gca, ['spots-' num2str(i) '.jpg']);
    end

    spots{i} = better_spots;
  end

  links = track_spots(spots, opts);
  [spots, links] = close_tracking_gaps(spots, links);

  centrosomes = get_struct('ruffles', [nframes, 1]);
  for i = 1:nframes
    centrosomes(i).carth = spots{i}(:, [2 1]);
    centrosomes(i).properties = spots{i}(:, 3:end);
    centrosomes(i).cluster = links{i};
  end

  mymovie.data.centrosomes = centrosomes;

  return;
end
