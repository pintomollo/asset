function mystruct = smooth_segmentation(mystruct, nimg, opts)

  if (nargin == 2)

    opts = nimg;

    if (isfield(mystruct, 'experiment'))
      error('Function not implemented yet...');
    else
      nframes = size_data(mystruct);
      for i=1:nframes
        mystruct = smooth_segmentation(mystruct, i, opts);
      end
    end

    return;
  end

  egg = mystruct.eggshell(nimg).carth;
  cortex = mystruct.cortex(nimg).carth;
  ruffles = mystruct.ruffles(nimg).carth;
  bounds = mystruct.ruffles(nimg).bounds;
  [pts, total_dist] = carth2linear(egg);
  pts = pts*total_dist;

  [imfs_x, imfs_y] = circular_emd(pts, egg, total_dist);
  egg = [imfs_x(:,end) imfs_y(:,end)];

  %plot(imfs_x(:,end), imfs_y(:,end), 'c')
  %plot(sum(imfs_x(:,end-1:end),2), sum(imfs_y(:,end-1:end),2), 'y')

  [pts, total_dist] = carth2linear(cortex);
  pts = pts*total_dist;
  [imfs_x, imfs_y] = circular_emd(pts, cortex, total_dist);

  if (size(imfs_x, 2) > 2)
    cortex = sum(imfs_x(:,end-1:end),2);
    sharper = sum(imfs_x(:,end-2:end),2);
  elseif (size(imfs_x, 2) == 2)
    cortex = imfs_x(:,end);
    sharper = sum(imfs_x(:,end-1:end),2);
  else
    cortex = imfs_x;
    sharper = imfs_x;
  end

  if (size(imfs_y, 2) > 2)
    cortex = [cortex, sum(imfs_y(:,end-1:end),2)];
    sharper = [sharper, sum(imfs_y(:,end-2:end),2)];
  elseif (size(imfs_y, 2) == 2)
    cortex = [cortex, imfs_y(:,end)];
    sharper = [sharper, sum(imfs_y(:,end-1:end),2)];
  else
    cortex = [cortex imfs_y];
    sharper = [sharper imfs_y];
  end


%  imshow(load_data(mymovie.cortex, nimg));
%  hold on
%  plot(egg(:,1), egg(:,2), 'g');
%  plot(cortex(:,1), cortex(:,2), 'm');
%  plot(sharper(:,1), sharper(:,2), 'k');
%  scatter(mymovie.markers.ruffles(nimg).carth(:,1), mymovie.markers.ruffles(nimg).carth(:,2), 'r')

%  new_cortex = cortex;

  for i=1:size(ruffles, 1)
    ruf = [ruffles(i,:); bounds(i,1:2); bounds(i,3:4)];
    dist = (bsxfun(@minus, ruf(:,1).', cortex(:,1)).^2 + bsxfun(@minus, ruf(:,2).', cortex(:,2)).^2);
    [junk, indx] = min(dist, [], 1);
    %scatter(cortex(indx,1), cortex(indx, 2),'b');
    
    width = pts(indx(3)) - pts(indx(2));
    if (width < 0)
      width = width + total_dist;
    end
    center = width/2 + pts(indx(2));

    circular_mask = fusion_mask(pts, center, 1.125*width/2, 0.1, total_dist);

    cortex = bsxfun(@times, sharper, circular_mask) + bsxfun(@times, cortex, (1-circular_mask));
  end

  mystruct.eggshell(nimg).carth = egg;;
  mystruct.cortex(nimg).carth = cortex;
  
%  plot(cortex(:,1), cortex(:,2), 'r');
  %drawnow;

  return;
end

function [imfs_x, imfs_y] = circular_emd(pos, pts, total_dist);

  
  circular_mask = fusion_mask(pos, total_dist/2, total_dist/3, 0.05, total_dist);
  imfs_x = emdc(pos, pts(:,1), true, 3).';
  imfs_y = emdc(pos, pts(:,2), true, 3).';

  half = round(length(pos)/2);
  pos = [pos(half+1:end, 1); pos(1:half,1)+total_dist];
  pts = [pts(half+1:end, :); pts(1:half,:)];

  tmp_x = emdc(pos, pts(:,1), true, 3).';
  tmp_y = emdc(pos, pts(:,2), true, 3).';

  tmp_x = [tmp_x(end-half+1:end,:); tmp_x(1:end-half,:)];
  tmp_y = [tmp_y(end-half+1:end,:); tmp_y(1:end-half,:)];

  [m] = size(imfs_x,2);
  [n] = size(tmp_x,2);

  if (m~=n)
    if (m>n)
      tmp_x = padarray(tmp_x, [0 m-n], 'pre');
    else
      imfs_x = padarray(imfs_x, [0 n-m], 'pre');
    end
  end

  [m] = size(imfs_y,2);
  [n] = size(tmp_y,2);

  if (m~=n)
    if (m>n)
      tmp_y = padarray(tmp_y, [0 m-n], 'pre');
    else
      imfs_y = padarray(imfs_y, [0 n-m], 'pre');
    end
  end

  imfs_x = bsxfun(@times, imfs_x, circular_mask) + bsxfun(@times, tmp_x, 1-circular_mask);
  imfs_y = bsxfun(@times, imfs_y, circular_mask) + bsxfun(@times, tmp_y, 1-circular_mask);

  return;
end

function [mask, width] = fusion_mask(pos, mu, sigma, slope, l)

  mask = zeros(size(pos));
  right = (pos > mu);

  slope = slope/2;
  while(any(mask(right) ~= 2))
    slope = slope*2;
    mask = (erf(slope*((pos-mu)+sigma)) + 1);
  end
  mask(right) = 0;

  if (mask(1)~=0)
    mask = mask + (erf(slope*((pos-mu-l)+sigma)) + 1);
  end
  
  mask(right) = mask(right) + (erf(-slope*((pos(right)-mu)-sigma)) + 1);
  if (mask(end)~=0)
    mask = mask + (erf(-slope*((pos-mu+l)-sigma)) + 1);
  end

  mask = mask*0.5;

  return;
end
