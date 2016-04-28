function new_img = imbin(img, bin_size)

  if (nargin == 1)
    bin_size = 2;
  end

  filter = ones(2*bin_size - 1);
  filter(1:bin_size-1, :) = 0;
  filter(:, 1:bin_size-1) = 0;
  filter = filter / sum(filter(:));

  img = imfilter(img, filter, 'symmetric');
  new_img = img(1:bin_size:end, 1:bin_size:end);

  return;
end
