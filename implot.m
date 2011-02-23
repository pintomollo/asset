function hhh = implot(img)

  if(~isnumeric(img))
    img = double(img);
  end

  full_img = zeros(size(img)+1);
  full_img(1:end-1,1:end-1) = img;

  [h,w] = size(img);
  [X,Y] = meshgrid([1:w+1]-0.5,[1:h+1]-0.5);

  hh=pcolor(X,Y,full_img);
  set(hh,'EdgeColor','none');
  %set(get(hh,'Parent'),'XAxisLocation','top','YDir','reverse');

  if(nargout>0)
    hhh=hh;
  end

  return;
end
