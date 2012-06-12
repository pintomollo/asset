function [new_path, shift] = detect_eggshell(img, path, axes_length, safety, pos, shift)

  if (nargin < 6)
    imgsize = size(img);

    edges = imadm_mex(img);
    inner = NaN(imgsize(1),1);

    for i=1:imgsize(1)
      tmp = find(edges(i,1:path(i,1)-1)>0, 1, 'last');
      if (~isempty(tmp))
        inner(i,1) = tmp;
      end
    end

    %figure;imshow(img);
    %hold on;
    %plot(path, 1:imgsize(1), 'b');
    %plot(inner, 1:imgsize(1), 'r');
  
    ell1 = pixels2elliptic(path, imgsize, axes_length, safety);
    ell2 = pixels2elliptic(inner, imgsize, axes_length, safety);
    ell0 = ell1;

    pts1 = elliptic2carth(ell1, [0;0], axes_length, 0);
    pts2 = elliptic2carth(ell2, [0;0], axes_length, 0);

    dist = sqrt(sum((pts1 - pts2).^2, 2));
    shift = median(dist(~isnan(dist)));

    norm = bsxfun(@rdivide, -pts1, sqrt(sum(pts1.^2, 2)));
    dpos = norm * pos * shift;

    %figure;plot([pts1(:,1) zeros(imgsize(1), 1)].', [pts1(:,2), zeros(imgsize(1), 1)].', 'k');
    %hold on;
    %scatter(pts1(:,1), pts1(:, 2), 'b');
    %scatter(pts2(:,1), pts2(:, 2), 'r');
    %scatter(pts1(:,1)+dpos(:,1), pts1(:, 2)+dpos(:,2), 'g');
    %scatter(pts1(:,1)+norm(:,1)*shift, pts1(:, 2)+norm(:,2)*shift, 'b');

    %ell1 = carth2elliptic(pts1, [0;0], axes_length, 0, 'radial');
    %correction = (1 - pos / sqrt((axes_length(1)*cos(ell1(:,1))).^2 + (axes_length(2)*sin(ell1(:,1))).^2));
    %ell1(:,2) = ell1(:,2) .* correction(:);
    %pts1 = elliptic2carth(ell1, [0;0], axes_length, 0, 'radial');

    ell1 = carth2elliptic(pts1 + dpos, [0;0], axes_length, 0);
    new_path = elliptic2pixels(ell1, imgsize, axes_length, safety);
  
    %plot(new_path(:,2), 1:imgsize(1), 'g');

    %figure;hold on;
    %myplot(ell0, 'b');
    %myplot(ell2, 'r');
    %myplot(ell1, 'g');
    %myplot(ell0(:,1), dist, 'k');
    %myplot(ell0(:,1), correction, 'm');

    %keyboard

%    theta = [0:2*pi/nrows:2*pi];
%    theta = theta(1:end-1).';
%
%    dist = path - inner;
%    correction = sqrt(cos(theta).^2 + (axes_length(2)/axes_length(1) * sin(theta)).^2);
%    pol_dist = dist .* correction;
%    best = median(pol_dist);
%
%    shift = best ./ correction;
%
%    new_path = path - pos*shift;
%
%    shift = best / axes_length(1);

  else

    if (prod(size(img)) > 2)
      imgsize = size(img);
    else
      imgsize = img;
    end
 
    ell1 = pixels2elliptic(path, imgsize, axes_length, safety);
    pts1 = elliptic2carth(ell1, [0;0], axes_length, 0);

    norm = bsxfun(@rdivide, -pts1, sqrt(sum(pts1.^2, 2)));
    dpos = norm * shift;

    %figure;plot([pts1(:,1) zeros(imgsize(1), 1)].', [pts1(:,2), zeros(imgsize(1), 1)].', 'k');
    %hold on;
    %scatter(pts1(:,1), pts1(:, 2), 'b');
    %scatter(pts2(:,1), pts2(:, 2), 'r');
    %scatter(pts1(:,1)+dpos(:,1)*(1-pos), pts1(:, 2)+dpos(:,2)*(1-pos), 'g');
    %scatter(pts1(:,1)-dpos(:,1)*(pos), pts1(:, 2)-dpos(:,2)*(pos), 'r');
    %scatter(pts1(:,1)+norm(:,1)*shift, pts1(:, 2)+norm(:,2)*shift, 'b');


    %ell1 = carth2elliptic(pts1, [0;0], axes_length, 0, 'radial');
    %correction = (ell1(:,2) ./ sqrt((axes_length(1)*cos(ell1(:,1))).^2 + (axes_length(2)*sin(ell1(:,1))).^2));

    %outer = ell1(:, 2) + pos*correction(:);
    %inner = ell1(:, 2) - (1-pos)*correction(:);

    %pts1 = elliptic2carth(ell1(:,1), inner, [0;0], axes_length, 0, 'radial');
    %pts2 = elliptic2carth(ell1(:,1), outer, [0;0], axes_length, 0, 'radial');

    ell1 = carth2elliptic(pts1 + dpos*(1-pos), [0;0], axes_length, 0);
    ell2 = carth2elliptic(pts1 - dpos*pos, [0;0], axes_length, 0);
    inner_path = elliptic2pixels(ell1, imgsize, axes_length, safety);
    outer_path = elliptic2pixels(ell2, imgsize, axes_length, safety);

   
    %theta = [0:2*pi/nrows:2*pi];
    %theta = theta(1:end-1).';

    %correction = sqrt(cos(theta).^2 + (axes_length(2)/axes_length(1) * sin(theta)).^2);

    %shift = shift * axes_length(1);
    %shift = shift ./ correction;

    %inner_path = path - (1 - pos) * shift;
    %outer_path = path + pos * shift;

    [new_path, shift] = deal(inner_path(:,2), outer_path(:,2));
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
