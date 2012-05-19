function [vars, ranges] = estimate_noise(img, wsize)

  if (nargin == 1)
    wsize = 5;
  end

  if (numel(unique(wsize)) > 1)

    vars = [];
    ranges = [];
    wsize = unique(wsize);

    for i=1:numel(wsize)
      [tmp_v, tmp_r] = estimate_noise(img, wsize(i));
      vars = [vars, tmp_v];
      ranges = [tmp_r, ranges];
    end

    return;
  end

  if (numel(wsize == 1))
    wsize = wsize([1 1]);
  end

  vars = [];
  
  ming = imfilter(img, fspecial('average', wsize), 'symmetric');
  lvar = (img - ming).^2;
  lvar(isnan(img)) = NaN;
  vars(end+1) = mymean(lvar(:));
  %ranges(1) = range(ming(:));
  
  ming = medfilt2(img, wsize);
  mvar = (img - ming).^2;
  mvar(isnan(img)) = NaN;
  vars(end+1) = mymean(mvar(:));
  %ranges(2) = range(ming(:));

  ming = imfilter(img, fspecial('gaussian', wsize, wsize(1)/3), 'symmetric');
  %ming = wiener2(img, wsize);
  gvar = (img - ming).^2;
  gvar(isnan(img)) = NaN;
  vars(end+1) = mymean(gvar(:));
  %ranges(3) = range(ming(:));

  ming = wiener2(img, wsize);
  wvar = (img - ming).^2;
  wvar(isnan(img)) = NaN;
  vars(end+1) = mymean(wvar(:));

  vars = sqrt(vars);
  ranges = range(img(:)) -8*vars;
  %ranges = [range(img(:)) - 8.5*vars; ...
  %          range(img(:)) - 7.5*vars];

  return;
end
