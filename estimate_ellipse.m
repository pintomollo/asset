function [center, axes_length, orientation, mask, estimation] = estimate_ellipse(img)

  'Use the paramters from opts instead'

  img = double(img);

  npixels = max(size(img));
  ntotal = prod(size(img));

  size250 = round(npixels/250);
  size200 = round(npixels/200);

  %img = imfilter(img,fspecial('average',round(npixels/100)),'symmetric');

  img = gaussian_mex(img, size250);
  img = median_mex(img, size200, 3);
  %for j=1:3
  %  img=medfilt2(img,[round(npixels/100) round(npixels/100)]);
  %end
  %figure;imshow(img)

  img = imnorm(img);
  %img = (img - min(min(img))) / (max(max(img)) -min(min(img)));

  %thresh = mean(mean(img));
  thresh = graythresh(img);
  img=(img>thresh);
  %figure;imshow(img)

  img = bwareaopen(img,size200);
  img = imclose(img, strel('disk', size250));
  %img = imdilate(img, ser5);
  %img = imfill(img,'holes');
  %img = imerode(img, ser5);
  img = bwareaopen(img,size200);

  img = separate(img);
  %figure;imshow(img)

  if (~any(any(img)))

    center = [];
    axes_length = [];
    orientation= [];
    mask = [];
    estimation = [];

    return;
  end

  [center, axes_length, orientation] = blob2ellipse(img);

  if (nargout > 3)
    [x,y] = draw_ellipse(center,axes_length,orientation);
    mask = roipoly(size(img,1),size(img,2),x,y);
  end

  if (nargout > 4)
    estimation = bwboundaries(img, 8, 'noholes');
    estimation = estimation{1};
    estimation = estimation(:,[2 1]);
  end

  %figure;imshow(img)
  %hold on;plot(x,y,'g', 'LineWidth',2)
 
  return;
end
