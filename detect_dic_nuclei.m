function mymovie =  detect_dic_nuclei(mymovie, opts)

  nuclei_size = 1;
  thresh = 0.75;
  area_thresh = 0.5;
  %thresh = [0.65:0.01:0.75];
  %morph_size = [1];

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

  %rand_frames = randperm(nframes);
  %classes = cell([10 1]);

  nuclei = get_struct('ruffles', [1, nframes]);

  for i=1:nframes
    nimg = i;
    %nimg = randi(nframes, 1);
    %nimg = 89;
    %nimg = rand_frames(i);

    %keyboard

    img = imnorm(double(load_data(mymovie.dic, nimg)));
    %imshow(img);
    %colormap('gray')
    %saveas(gca, ['frame-' num2str(nimg) '-ref.png']);
    %colormap('jet')

    %for n=nuclei_size
    %    for m=morph_size


    entr = imnorm(entropyfilt(img, filt));

    path = mymovie.dic.cortex(nimg).carth;
    mask = roipoly(blank_img, path(:,1), path(:,2));

    %entr(~mask) = NaN;
    %entr = imnorm(entr);
    %thresh = graythresh(entr) * t;
    entr(~mask) = 1;

    %  for t=thresh
    bw = (entr < thresh);
    %bw(~mask) = 0;
    %bw = imclose(bw, strel('disk', ceil(filt_size/m)));
    %stats = regionprops(bw, img, 'Area', 'Eccentricity', 'MeanIntensity', 'MajorAxisLength');

    orig_bw = bw;
    
    %for a=area_thresh

    areas = round(pi * (area_thresh / opts.pixel_size)^2);

    bw = bwareaopen(orig_bw, areas);
    bw = imopen(bw, strel('disk', ceil(filt_size)));
    bw = imfill(bw, 'holes');

    stats = regionprops(bw, 'MajorAxisLength', 'Centroid', 'Area');

      %imshow(img);
      %hold on
    nucleus = get_struct('ruffles');
    for r=1:length(stats)
      if (stats(r).MajorAxisLength >= 50)
      %  continue;
      %  worse = [worse; [stats(r).Centroid, stats(r).Area, nimg]];
      %else
      
      nucleus.carth = [nucleus.carth; stats(r).Centroid];
      nucleus.properties = [nucleus.properties; stats(r).Area];

        %goods = [goods; [stats(r).Centroid, stats(r).Area, nimg]];
      end

      %box_pos = stats(r).BoundingBox(1, 1:2);
      %box_width = stats(r).BoundingBox(1, 3:4);

      %plot(box_pos(1)+[0 0; box_width([1 1]); 0 box_width(1); 0 box_width(1)], [0 box_width(2); 0 box_width(2); 0 0; box_width([2 2])]+box_pos(2), 'k')



      %spot_type = input('') + 1;
      %classes{spot_type} = [classes{spot_type}; [stats(r).Eccentricity, stats(r).MeanIntensity, stats(r).MajorAxisLength, stats(r).Solidity, stats(r).Area]];
    end
      %hold off;
      %pause

    %save('spots.mat', 'classes');

    l = bwlabel(bw, 8);
    l(l==0) = img(l==0);
    imagesc(l);
    %drawnow
    %subplot(121);imshow(img);
    %subplot(122);imagesc(bwlabel(bw, 8));
    title(num2str(nimg));

    saveas(gca, ['frame-' num2str(nimg) '.png']);

    %end
%    end
%    end
    nuclei(nimg) = nucleus;
  end

  %mymovie.dic.nuclei = nuclei;
  mymovie.dic.nuclei = tracking(nuclei, opts);

  return;
end
