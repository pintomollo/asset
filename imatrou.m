function atrous = imatrou(img, size_max)

  [h,w] = size(img);
  if (nargin < 2)
    size_max = max(h,w);
  end
  kernel = [1 4 6 4 1] / 16;

  nplanes = floor(log2(size_max - 1) - 1);
  atrous = zeros(h,w,nplanes+1);

  for i=1:nplanes
    atrou_kernel = insert_holes(kernel, i);
    filtered_img = imfilter(img, atrou_kernel, 'symmetric');
    filtered_img = imfilter(filtered_img, atrou_kernel.', 'symmetric');

    atrous(:,:,i) = img - filtered_img;
    img = filtered_img;

    %figure;implot(img);
  end
  atrous(:,:,end) = img;

  return;
end

function atrou_kernel = insert_holes(kernel, i)

  if (i <= 1)
    atrou_kernel = kernel;

    return;
  end

  ntaps = length(kernel);
  ninserts = (2^(i-1) - 1);
  new_size = (ntaps-1)*ninserts + ntaps;

  atrou_kernel = zeros(1, new_size);
  atrou_kernel(1:ninserts+1:end) = kernel;

  return;
end
