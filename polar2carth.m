function [carth_pts] = polar2carth(pts, imgsize, center)

  [s1 s2] = size(pts);

  if(s1==1)
    pts = pts';
    [s1 s2] = size(pts);
  end
  if(size(center,2)==1)
    center=center';
  end

  if(s2==1)
    O = [0:2*pi/length(pts):2*pi]';
    O = O(1:end-1);
    s2=2;
  else
    O = pts(:,1);
    if(max(O)>2*pi)
      O = O * (2*pi/imgsize(1));
    end
  end

  r = pts(:,2) * (sqrt(imgsize(2)^2+(imgsize(1)/2)^2)/(2*imgsize(2)));

  carth_pts = zeros([s1 s2]);
  carth_pts(:,1) = r.*cos(O) + center(1);
  carth_pts(:,2) = r.*sin(O) + center(2);

  return;
end
