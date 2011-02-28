function img = imhotpixels(orig_img, thresh, method)

  if (nargin == 1)
    thresh = 15;
    method = 'median';
  elseif (nargin == 2)
    if (ischar(thresh))
      method = thresh;
      thresh = 15;
    else
      method = 'median';
    end
  end

  pixels = double(orig_img(:));
  mean_value = mean(pixels);
  stddev = std(pixels);

  thresh = thresh * stddev;

  bad_pixels = ((orig_img > (mean_value + thresh)) | (orig_img < mean_value - thresh));

  img = orig_img;

  if (any(any(bad_pixels)))
  
    switch (method)
      case 'median'
        %filt_img = medfilt2(orig_img);
        filt_img = median_mex(double(orig_img));
      otherwise
        filt_img = imfilter(double(orig_img), fspecial(method), 'symmetric');
    end

    img(bad_pixels) = filt_img(bad_pixels);
  end

end
