function [all_estim] = detect_embryos(imgs, use_edges, egg_min_size, intense_thresh, border_thresh)

  [h,w,nframes] = size(imgs);
  imgsize = [h, w];

  rad_thresh = max(egg_min_size);
  rad_thresh = rad_thresh ./ [0.5 4 8 10 12];
  area_thresh = prod(egg_min_size + rad_thresh(3))*pi;

  all_estim = cell(nframes, 1);

  for nming = 1:nframes
    img = double(imgs(:,:,n));

    if (use_edges)
      img = imadm_mex(img);
      thresh = graythresh(img);
      img = (img > thresh*intense_thresh*(max(img(:))));
    else
      img = gaussian_mex(img, rad_thresh(5));
      img = median_mex(img, rad_thresh(4), 3);

      img = (img > intense_thresh);
    end

    img = padarray(img, rad_thresh([1 1]));

    img = imdilate(img, strel('disk', rad_thresh(3)));
    img = imfill(img, 'holes');
    img = bwareaopen(img, area_thresh);
    img = imdilate(img, strel('disk', rad_thresh(2)));
    img = imfill(img, 'holes');
    img = imerode(img, strel('disk', rad_thresh(2)+rad_thresh(3)));

    img =  img([1:imgsize(1)]+rad_thresh(1),[1:imgsize(2)]+rad_thresh(1));

    img =  img((size10+1):(size10+imgsize(1)),(size10+1):(size10+imgsize(2)));

    if (opts.verbosity == 3)
      figure;imshow(img);
    end

    if(~any(img) | (sum(img(:)) / prod(imgsize) > 0.9))
      continue;
    end

    estim = bwboundaries(img, 8, 'noholes');

    if (isempty(estim))
      continue;
    end

    estim = cellfun(@(x)([x; NaN(1,2)]), estim, 'UniformOutput', false);
    estim = cat(1, estim{:});

    borders = (any(estim <= border_thresh | bsxfun(@eq, estim, imgsize([2 1])-border_thresh+1), 2));

    all_estim{nimg} = estim(~borders, :);
  end

  if (nframes == 1)
    all_estim = all_estim{1};
  end

  return;
end
