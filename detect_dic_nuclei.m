function mymovie =  detect_dic_nuclei(mymovie, opts)

  nuclei_size = 1;
  thresh = 0.75;
  area_thresh = 0.5;

  filt_size = (nuclei_size / opts.pixel_size);
  if (mod(filt_size, 2) < 1)
    filt_size = ceil(filt_size);
  else
    filt_size = floor(filt_size);
  end
  filt = true(filt_size);

  [nframes, imgsize] = size_data(mymovie.dic);
  blank_img = false(imgsize);

  goods = zeros(0, 4);
  worse = zeros(0, 4);

  nuclei = get_struct('ruffles', [nframes, 1]);

  for i=1:nframes
    nimg = i;

    img = imnorm(double(load_data(mymovie.dic, nimg)));
    entr = imnorm(entropyfilt(img, filt));

    path = mymovie.dic.cortex(nimg).carth;
    mask = roipoly(blank_img, path(:,1), path(:,2));
    entr(~mask) = 1;

    bw = (entr < thresh);
    orig_bw = bw;

    areas = round(pi * (area_thresh / opts.pixel_size)^2);

    bw = bwareaopen(orig_bw, areas);
    bw = imopen(bw, strel('disk', ceil(filt_size)));
    bw = imfill(bw, 'holes');

    stats = regionprops(bw, 'MajorAxisLength', 'Centroid', 'Area');
    nucleus = get_struct('ruffles');
    for r=1:length(stats)
      if (stats(r).MajorAxisLength >= 50)
        if (isnan(nucleus.carth(1)))
          nucleus.carth = stats(r).Centroid;
        else
          nucleus.carth = [nucleus.carth; stats(r).Centroid];
        end
        nucleus.properties = [nucleus.properties; stats(r).Area];
      end
    end

    %l = bwlabel(bw, 8);
    %l(l==0) = img(l==0);
    %imagesc(l);
    %title(num2str(nimg));

    %saveas(gca, ['frame-' num2str(nimg) '.png']);
    nuclei(nimg) = nucleus;
  end

  mymovie.dic.nuclei = tracking(nuclei, opts);

  return;
end
