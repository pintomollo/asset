function [new_img] = imcosmics(img, block_size, thresh)
% IMCOSMICS implements the cosmic rays removal from [1]. This algorithm works
% by removing pixels which values are separated from lower values by a gap larger
% than a threshold estimated using the median and the mad.
%
%   IMG = IMCOSMICS(RAYS, BLOCK_SIZE, THRESH) removes the cosmic rays in RAYS and
%   returns the cleaned IMG. The algorithm works by processing contiguous tiles of
%   size [BLOCK_SIZE BLOCK_SIZE] with an additional overlap of 0.5*BLOCK_SIZE in both
%   dimensions. Cosmic ray pixels are identified in each tile as the pixels which
%   value are separated from lower values by a gap larger than THRESH*MAD. These
%   pixels are then replaced by the median value in the corresponding tile. This
%   process is repeated iteratively until convergeance of the detection (up to 5 times).
%
%   IMG = IMCOSMICS(RAYS, BLOCK_SIZE) utilizes the published THRESH value of 3.
%
%   IMG = IMCOSMICS(RAYS) utilizes in addition the published BLOCK_SIZE value of 10.
%
%   IMGS = IMCOSMICS(STACK, ...) processes each plane of the stack independently.
%
% References:
%   [1] Pych, W. A Fast Algorithm for Cosmic-Ray Removal from Single Images.
%       Publications of the Astronomical Society of the Pacific 816 (2004).
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 15.05.2014

  % The default values in case of insufficient inputs
  if (nargin == 1)
    block_size = 10;
    thresh = 3;
  elseif (nargin == 2)
    thresh = 3;
  end

  % No data to analyze
  if (nargin == 0 || isempty(img))
    new_img = [];
    return;
  end

  % In case of a stack, simply process each slice independently
  nplanes = size(img, 3);
  if (nplanes > 1)

    % Allocate the output and iterate
    new_img = img;
    for i = 1:nplanes
      new_img(:,:,i) = imcosmics(img(:,:,i), block_size, thresh);
    end

    return;
  end

  % Get the type of img as we need to convert it to double for the division
  class_type = class(img);
  is_double = (class_type(1) == 'd');

  % In that case, convert it
  if (~is_double)
    img = double(img);
  end

  % Refresh the image package
  %pkg unload image;
  %pkg load image;
  %pkg unload statistics;
  %pkg load statistics;

  % The amount of overlap and the indexes to extract the center of the block
  overlap = ceil(block_size/2);
  bindex = [1:block_size] + overlap;
  iindex = [1:size(img,1)];
  jindex = [1:size(img,2)];

  % In case of a single plane, [1] advises to process iteratively the filtering
  % up to 5 times, even though convergeance is usually reached in 2-3 calls.
  for i = 1:5

    % We call blockproc to obtain the block partitioning of the image, adding to that
    % an overlap of block_size/2. The filtering is done by filter_cosmics
    new_img = blockproc(img, [block_size block_size], [overlap overlap], @filter_cosmics);
    new_img = new_img(iindex, jindex);

    % Check for convergeance
    if (all(new_img(:) == img(:)))
      break;
    end

    % Otherwise, loop
    img = new_img;
  end

  % And convert back to the original type
  if (~is_double)
    new_img = cast(new_img, class_type);
  end

  %pkg unload image;
  %pkg unload statistics;

  return;

  function res = filter_cosmics(block_struct)
  % This nested function performs the actual filtering, having access to all the
  % required parameters directly from the main function.

    % Extract the data from the blockproc structure.
    res = block_struct;

    % We do not need to have them ordered spatially
    data = res(:);

    % 0 padding indicates pixels outside from the image. We assume that no real
    % pixel will have this exact value.
    data(data == 0) = NaN;

    % We need to robustly estimate the mean and the standard deviation of the pixel
    % distribution. In opposition with what is proposed in [1], we utilize the median
    % and the MAD estimators for this purpose.

    dmed = nanmedian(data);
    dmad = 1.4826 * nanmedian(abs(data-dmed));

    % Here we most likely have a problem, so try another approach to estimate the
    % standard deviation.
    if (dmad==0)
      dmad = 1.4826 * nanmean(abs(data-dmed));

      % If nothing else worked, follow [1] to esimate the standard deviation
      if (dmad == 0)
        dmad = sqrt(nanmean(data.^2 - dmed^2));
      end
    end

    % Sort the data, keeping track of the correspondency
    [data, indx] = sort(data);

    % And measure the distance between pixel values
    dist = [0; diff(data)];

    % Our cosmic rays are higher than the estimated mean and have a gap larger than
    % the standard deviation times the threshold
    bads = (data > dmed & dist > thresh*dmad);

    % If we find one, [1] removes all the higher-valued pixels. Given that we have
    % ordere them, we find the first index and replace all higher values by the median
    if (any(bads))
      first = find(bads, 1, 'first');
      res(indx(first:end)) = dmed;
    end

    res = res(bindex, bindex);

    return;
  end
end
