function path = adapt_path(imgsize, pts)

  if (prod(size(imgsize)) > 2)
    imgsize = size(imgsize);
  end

  if (size(pts,2)==1)
    pos = [1:imgsize(1)/(length(pts)-1):imgsize(1)];
    pts = [pos.' pts];
    %theta = [0:2*pi/length(pts):2*pi];
    %pts = [theta(1:end-1)' pts];
  else
    [junk, indx] = sort(pts(:,1));
    pts = pts(indx,:);
  end

  [npts, junk] = size(pts);
  full_pos = [1:imgsize(1)].';
  %full_theta = [0:2*pi/imgsize(1):2*pi];
  %full_theta = full_theta(1:end-1)';

  start = find(pts(:,1) < full_pos(end) - full_pos(1), 1, 'last');
  ends = find(pts(:,1) > full_pos(1), 1, 'first');
  
  pts = pts([start:end 1:end 1:ends], :);
  pts(1:npts-start+1,1) = pts(1:npts-start+1,1) - full_pos(end);
  pts(end-ends+1:end,1) = pts(end-ends+1:end,1) + full_pos(end);

  path = interp1q(pts(:,1), pts(:,2), full_pos);
  %tmpegg = pts;
  %tmpegg = tmpegg([end 1:end 1],:);
  %tmpegg(1,1) = tmpegg(1,1) - 2*pi;
  %tmpegg(end,1) = tmpegg(end,1) + 2*pi;

  %path = interp1q(tmpegg(:,1),tmpegg(:,2),full_theta);
  %path = path * imgsize(2) / safety;

  return;
end
