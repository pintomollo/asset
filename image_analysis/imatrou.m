function [projection, atrous] = imatrou(img, size_max, coef)
% IMATROU implements the spot detection algorithm from [1]. This algorithm
%   works by recursively filtering the image with kernels of increasing size. This
%   has the advantage of detecting objects with any size.
%
%   SPOTS = IMATROU(IMG, SPOT_SIZE, NOISE_THRESH) performs the "a trous" wavelet
%   decomposition of IMG as described in [1]. It will decompose IMG with separable
%   kernels up to SPOT_SIZE. Note that the first kernel is [1/16 1/4 3/8 1/4 1/16].
%   In addition, noise removal is performed as described in [1] by removing hits 
%   that are smaller than NOISE_THRESH * MAD(IMG) / 0.67. SPOTS is then defined as
%   the multi-scale product of this filtered decomposition as actual spots are 
%   correlated across resolution levels, while noise is not. Note that setting 
%   NOISE_THRESH to 0 will remove the noise filtering. Note also that planes that
%   have no detected spots are set to 1 in order not tocorruptt the product.
%
%   SPOTS = IMATROU(IMG, SPOT_SIZE) uses the default value of 3 for NOISE_THRESH 
%   as proposed in [1].
%
%   SPOTS = IMATROU(IMG) uses the default value MAX(SIZE(IMG)) for SPOT_SIZE.
%
%   [SPOTS, DECOMPOSITION] = IMATROU(...) returns in addition the wavelet 
%   decomposition as a 3D matrix. Note that the last plane of the stack is the
%   final smoothed image which should not be used for the multi-scale product.
%
% References:
%   [1] Olivo-Marin, J.-C. Extraction of spots in biological images using 
%       multiscale products. Pattern Recognit. 35, 1989-1996 (2002).
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 13.01.2011

  % Get the size of the image
  [h, w] = size(img);

  % Check the input arguments
  if (nargin < 2)
    coef = 3;
    size_max = max(h, w);
  elseif (nargin < 3)
    coef = 3;
  end

  % Create the kernel
  kernel = [1 4 6 4 1];
  kernel_norm = 0.0625^2; %1/16;

  % Compute the number of planes that we'll create
  nplanes = floor(log2(size_max - 1) - 1);

  % Initialize the decomposition
  atrous = ones(h,w,nplanes+1);

  if (coef > 0)
    % Loop over the different sizes of kernel
    for i = 1:nplanes

      % Create the kernel of the correct size
      atrou_kernel = insert_holes(kernel, i);

      % Filter the image sequentially in both dimension
      filtered_img = imfilter(img, atrou_kernel.', 'replicate');
      filtered_img = kernel_norm * imfilter(filtered_img, atrou_kernel, 'replicate');

      % The new plan (termed "deltail image") is defined as the difference 
      % between two consecutive levels of filtering (termed "smoothed approximation").
      decomposition = img - filtered_img;

      % Noise as defined in [1] using the MAD estimate of the detail plane
      noise = (abs(decomposition) < coef * 1.4826 * mad(decomposition(:), 1));

      % No information would kill the projection, consequently we keep the ones
      if (any(~noise(:)))

        % Otherwise we remove the noise and store the plane
        decomposition(noise) = 0;
        atrous(:,:,i) = decomposition;
      end

      % We reassign the filtered image to iteratively filter it
      img = filtered_img;
    end

    img = sum(atrous(:,:,1:end-1), 3);
  end

  % Loop over the different sizes of kernel
  for i = 1:nplanes

    % Create the kernel of the correct size
    atrou_kernel = insert_holes(kernel, i);

    % Filter the image sequentially in both dimension
    filtered_img = imfilter(img, atrou_kernel.', 'replicate');
    filtered_img = kernel_norm * imfilter(filtered_img, atrou_kernel, 'replicate');

    % The new plan (termed "detail image") is defined as the difference 
    % between two consecutive levels of filtering (termed "smoothed approximation").
    decomposition = img - filtered_img;

    % Otherwise we remove the noise and store the plane
    atrous(:,:,i) = decomposition;

    % We reassign the filtered image to iteratively filter it
    img = filtered_img;
  end

  % Compute the multi-scale product
  projection = abs(prod(atrous, 3));

  % If we need to return the decomposition, we store the final aproximation too
  if (nargout > 1)
    atrous(:,:,end) = img;
  else
    atrous = [];
  end

  return;
end

% This function creates the kernel of level i
function atrou_kernel = insert_holes(kernel, i)

  % Initially we keep the original kernel
  if (i <= 1)
    atrou_kernel = kernel;

    return;
  end

  % Count the number of elements ("taps") of the kernel
  ntaps = length(kernel);

  % And we insert that many zeros in between each of them (see [1])
  ninserts = (2^(i-1) - 1);

  % Final size of the new kernel
  new_size = (ntaps-1)*ninserts + ntaps;

  % Initiallize the kernel
  atrou_kernel = zeros(1, new_size);

  % Set the taps at the correct positions
  atrou_kernel(1:ninserts+1:end) = kernel;

  return;
end
