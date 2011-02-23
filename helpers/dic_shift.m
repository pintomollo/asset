function [best_coefs, shift, width] = dic_shift(mymovie, trackings)
%  Orig. Err ? std, Corr. Err ? std (% of Embryo radius)
%      0.0424    0.0312    0.0193    0.0174
%      a =
%          0.0424
%          -0.0440
%          -0.0176
%      b =
%
%          5.1836

  close all

  [imgsize, nframes] = size_data(mymovie.dic);

  centers = trackings.dic.reference.center;
  axes_length = trackings.dic.reference.axes_length;
  orientations = trackings.dic.reference.orientation;

  dic = trackings.dic.mean;
  markers = trackings.markers.mean;

  safety = 1.2;
  range = 0;
  corr = zeros(nframes,1);
  thetas = [0:pi/100:2*pi]';
  thetas = thetas(1:end-1);

  dic_path = zeros(length(thetas), nframes, 2);
  fluo_path = zeros(length(thetas), nframes, 2);
  intensity = zeros(length(thetas), nframes, range);

  for nimg = 1:nframes
    img = load_data(mymovie.dic,nimg);
    polar_img = elliptic_coordinate(img, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), safety);

    egg_path = carth2elliptic(dic(1,nimg), centers(:,nimg), axes_length(:,nimg), orientations(1,nimg),true);
    cortex_path = carth2elliptic(dic(2,nimg), centers(:,nimg), axes_length(:,nimg), orientations(1,nimg),true);
    fegg_path = carth2elliptic(markers(1,nimg), centers(:,nimg), axes_length(:,nimg), orientations(1,nimg),true);
    fcortex_path = carth2elliptic(markers(2,nimg), centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), true);

    tmp_path = fnval(egg_path, egg_path.breaks);
    egg_path = [egg_path.breaks' tmp_path(1,:)'];
    tmp_path = fnval(cortex_path, cortex_path.breaks);
    cortex_path = [cortex_path.breaks' tmp_path(1,:)'];
    tmp_path = fnval(fegg_path, fegg_path.breaks);
    fegg_path = [fegg_path.breaks' tmp_path(1,:)'];
    tmp_path = fnval(fcortex_path, fcortex_path.breaks);
    fcortex_path = [fcortex_path.breaks' tmp_path(1,:)'];

    pix_path = elliptic2pixels(egg_path, size(polar_img), safety);
    full_pixs = [1:size(polar_img,1)]';
    pix_path = adapt_path(size(polar_img), pix_path);

    %int = zeros(size(full_pixs));
    full_thetas = (full_pixs - 1) * 2*pi/full_pixs(end);
    for d=-range:range
      int = bilinear_mex(polar_img, pix_path+d, full_pixs);
      [int] = interp_path(thetas, [full_thetas, int]);
      intensity(:,nimg,d+range+1) = int(:,2);
    end
    %int = int / length(range);

    %dist = abs(full_thetas - x(3));
    %shift_indx = find(dist == min(dist));
    %int = int([shift_indx:end 1:shift_indx-1]);

    %[int] = interp_path(thetas, [full_thetas, int]);
    [egg_path] = interp_path(thetas, egg_path);
    [cortex_path] = interp_path(thetas, cortex_path);
    [fegg_path] = interp_path(thetas, fegg_path);
    [fcortex_path] = interp_path(thetas, fcortex_path);

    dic_path(:,nimg,1) = egg_path(:,2);
    fluo_path(:,nimg,1) = fegg_path(:,2);
    dic_path(:,nimg,2) = cortex_path(:,2);
    fluo_path(:,nimg,2) = fcortex_path(:,2);
  end

  dic_path = dic_path(:,:,2);
  fluo_path = fluo_path(:,:,2);

  min_se = Inf;
  shift = 0;
  width = 0;
  best_coefs = [];
  best_stats = [];

  corrs = dic_path - fluo_path;


  tmp = zeros(length(thetas), range+1);
  for j = 0:range
    int = mean(intensity(:,:,range+1+[-j:j]),3);
    %avgs = mean(int);
    mins = min(int);
    maxs = max(int)-mins;
    int = bsxfun(@minus, int, mins);
    int = bsxfun(@rdivide, int, maxs);
    %int = repmat(int, [1 1 2]);
    int = int(:);
    %avgs = repmat(avgs, [length(thetas) 1 2]);
    %mins = repmat(mins, [length(thetas) 1 2]);
    maxs = repmat(maxs, [length(thetas) 1 1]);
    %int = [int mins(:) maxs(:)];
    int = [int maxs(:)];

    for i=1:length(thetas)
      [coefs, stats] = robustfit(int, corrs(:));
      if (stats.s < min_se)
        min_se = stats.s;
        shift = i;
        best_coefs = coefs;
        best_stats = stats;
        width = j;
      end
      tmp(i,j+1) = stats.s;
      corrs = corrs([2:end 1],:,:);
    end
  end

  if (range > 0)
    [X,Y] = meshgrid([0:range],thetas);
    figure;surf(X,Y,tmp);
    hold on;
    scatter3(width, thetas(shift), min_se, 'r');
    xlim([0 range]);
    ylim([0 2*pi]);
    set(gca, 'Box', 'on', 'YTick', [0:pi/2:2*pi]);
  else
    figure;plot(thetas,tmp);
    hold on;
    scatter(thetas(shift), min_se, 'r');
    xlim([0 2*pi]);
    set(gca, 'Box', 'on', 'XTick', [0:pi/2:2*pi]);
  end

  int = mean(intensity(:,:,range+1+[-width:width]),3);
  %avgs = mean(int);
  mins = min(int);
  maxs = max(int)-mins;
  int = bsxfun(@minus, int, mins);
  int = bsxfun(@rdivide, int, maxs);
  %int = repmat(int, [1 1 2]);
  %avgs = repmat(avgs, [length(thetas) 1 2]);
  %mins = repmat(mins, [length(thetas) 1 2]);

  figure;
  [hax] = plotyy(thetas, cat(2, corrs(:,1:15), mean(corrs(:,1:15),2)), thetas, cat(2, int(:,1:15,:), mean(int(:,1:15,:),2)));
  set(hax, 'XLim', [0 2*pi], 'Box', 'on', 'XTick', [0:pi/2:2*pi]);


  corrs = corrs([shift:end 1:shift-1],:,:); 
  maxs = repmat(maxs, [length(thetas) 1 1]);
  init_dist = path_error(thetas, corrs, 36);

  y = corrs;
  f = y - (int*best_coefs(2) + best_coefs(1) + best_coefs(3)*maxs);
  y = y - mean(y(:));

  R2 = 1 - sum(f(:).^2)/sum(y(:).^2)

  corrs = corrs - best_coefs(3)*maxs;
  %int = int(:);;
  

  %figure;hold on;
  %scatter(int(:),corrs(:),'k');
  %plot([0 1],[0 1]*best_coefs(2) + best_coefs(1), 'r');
  %plot([0 1],[0 1]*(best_coefs(2)+best_stats.se(2)) + best_coefs(1) + best_stats.se(1), 'r');
  %plot([0 1],[0 1]*(best_coefs(2)-best_stats.se(2)) + best_coefs(1) - best_stats.se(1), 'r');
  %set(gca, 'Box', 'on')


  int = int([end-shift+2:end 1:end-shift+1],:,:);
  %avgs = avgs([end-shift+2:end 1:end-shift+1],:,:);
  %mins = mins([end-shift+2:end 1:end-shift+1],:,:);
  %maxs = maxs([end-shift+2:end 1:end-shift+1],:,:);
  expect = dic_path - (int*best_coefs(2) + best_coefs(1) + best_coefs(3)*maxs);
  improved_dist = path_error(thetas, expect - fluo_path, 36);

  [ferr, fstds] = mymean(init_dist(:));
  [ierr, istds] = mymean(improved_dist(:));
  
  disp('  Orig. Err ± std, Corr. Err ± std (% of Embryo radius)');
  disp([ferr, fstds, ierr, istds] * 100)

  %keyboard

  nimg=32;                                                                      
  img = load_data(mymovie.dic, nimg);                                           
  %cart_egg = elliptic2carth(thetas, expect(:,nimg, 1),centers(:,nimg),axes_length(:,nimg),orientations(1,nimg));
  cart_cortex = elliptic2carth(thetas, expect(:,nimg),centers(:,nimg),axes_length(:,nimg),orientations(1,nimg));

  rescale_size = [388 591];
  figure;hold on;
  imshow(realign(img,rescale_size,centers(:,nimg),orientations(1,nimg)));
  %myplot(realign(dic(1,nimg),rescale_size,centers(:,nimg),orientations(1,nimg)));
  myplot(realign(dic(2,nimg),rescale_size,centers(:,nimg),orientations(1,nimg)));        
  %myplot(realign(markers(1,nimg),rescale_size,centers(:,nimg),orientations(1,nimg)),'r');
  myplot(realign(markers(2,nimg),rescale_size,centers(:,nimg),orientations(1,nimg)),'r');
  %myplot(realign(cart_egg,rescale_size,centers(:,nimg),orientations(1,nimg)),'g');
  myplot(realign(cart_cortex,rescale_size,centers(:,nimg),orientations(1,nimg)),'g');

  tmp = init_dist(:,nimg,:);
  init_error = mymean(tmp(:)) * 100
  tmp = improved_dist(:,nimg,:);
  corr_error = mymean(tmp(:)) * 100

  %myplot(realign(mymovie.dic.cortex(nimg).carth,rescale_size,centers(:,nimg),orientations(1,nimg)),'m');
  %myplot(realign(mymovie.markers.cortex(nimg).carth,rescale_size,centers(:,nimg),orientations(1,nimg)),'y');

%  img = load_data(mymovie.dic,nimg);
%  polar_img = elliptic_coordinate(img, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), safety);
%
%  egg_path = carth2elliptic(mymovie.dic.eggshell(nimg).carth, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg));
%  cortex_path = carth2elliptic(mymovie.dic.cortex(nimg).carth, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg));
%  [cortex_path] = interp_path(thetas, cortex_path);
%  cortex_path = cortex_path(:,2);
%
%  pix_path = elliptic2pixels(egg_path, size(polar_img), safety);
%  pix_path = adapt_path(size(polar_img), pix_path);
%
%  full_pixs = [1:size(polar_img,1)]';
%  full_thetas = (full_pixs - 1) * 2*pi/full_pixs(end);
%  int = bilinear_mex(polar_img, pix_path, full_pixs);
%  [int] = interp_path(thetas, [full_thetas, int]);
%  int = int([end-shift+2:end 1:end-shift+1],2);
%
%  maxs = max(int);
%  expect = cortex_path(:) - (int(:)*best_coefs(2) + best_coefs(1) + best_coefs(3)*maxs);
%  expect = elliptic2carth(thetas, expect,centers(:,nimg),axes_length(:,nimg),orientations(1,nimg));
%  myplot(realign(expect,rescale_size,centers(:,nimg),orientations(1,nimg)),'w');
%
%  keyboard

  shift = thetas(shift);

  return;
end

function [new_ref, new_path] = interp_path(ref, path)

  if (nargout == 1)
    theta = ref(:,1);
  else
    theta = unique([ref(:,1); path(:,1)]);
  end

  if (all(path(end,:)==path(1,:)))
    path = path(1:end-1,:);
  end
  
  path = path([end 1:end 1],:);
  path(1,1) = path(1,1) - 2*pi;
  path(end,1) = path(end,1) + 2*pi;

  path = interp1q(path(:,1), path(:,2), theta);
  new_ref = [theta path];

  if (nargout == 1)
    return;
  end
  new_path = new_ref;

  if (all(ref(end,:)==ref(1,:)))
    ref = ref(1:end-1,:);
  end
  
  ref = ref([end 1:end 1],:);
  ref(1,1) = ref(1,1) - 2*pi;
  ref(end,1) = ref(end,1) + 2*pi;

  ref = interp1q(ref(:,1), ref(:,2), theta);
  new_ref = [theta ref];

  return;
end
