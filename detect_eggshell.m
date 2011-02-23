function [new_path, shift] = detect_egghell(img, path, axes_length, pos, shift)

  if (nargin < 5)
    [nrows,npts] = size(img);

    edges = imadm_mex(img);
    inner = zeros(nrows,1);

    for i=1:nrows
      tmp = find(edges(i,1:path(i,1)-1)>0, 1, 'last');
      if (~isempty(tmp))
        inner(i,1) = tmp;
      end
    end

    theta = [0:2*pi/nrows:2*pi];
    theta = theta(1:end-1).';

    dist = path - inner;
    correction = sqrt(cos(theta).^2 + (axes_length(2)/axes_length(1) * sin(theta)).^2);
    pol_dist = dist .* correction;
    best = median(pol_dist);

    shift = best ./ correction;
    %inner_path = path - shift;

    new_path = path - pos*shift;

    shift = best / axes_length(1);

  else

    if (prod(size(img)) > 2)
      [nrows, npts] = size(img);
    else
      nrows = img(1);
      npts = img(2);
    end
    
    theta = [0:2*pi/nrows:2*pi];
    theta = theta(1:end-1).';

    correction = sqrt(cos(theta).^2 + (axes_length(2)/axes_length(1) * sin(theta)).^2);

    shift = shift * axes_length(1);
    shift = shift ./ correction;

    inner_path = path - (1 - pos) * shift;
    outer_path = path + pos * shift;

    [new_path, shift] = deal(inner_path, outer_path);
  end

  
  %shift = median(pol_dist) ./ correction;
  %median_inner = path - shift;

  %plot(inner,[1:length(path)],'g')
  %plot(median_inner,[1:length(path)],'r')
  %plot(median_inner2,[1:length(path)],'m')

  %figure;imshow(new_img)
  %hold on;
  %plot(path,[1:length(path)])
  %plot(starts,[1:length(path)],'g')

  %figure;imshow(edges)
  %figure;imshow(imadm(new_img))

  %figure;hist(dist)
  %figure;hist(pol_dist,100)
  %kk

  return;
end
