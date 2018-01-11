function plot_3d_ellipse(pts, axes, orient)
  hold on;

  if (nargin == 3)
    centers = pts;

    nslices = 7;
    z_pos = [-nslices:nslices] * (axes(3)) / (nslices+1) ;
    z_coef = sqrt(1 - (z_pos.^2 / axes(3).^2));
    z_pos = z_pos + centers(3);

    for i=1:length(z_pos)
      [x, y] = draw_ellipse(centers(1:2), axes(1:2)*z_coef(i), orient);
      plot3(x,y,z_pos(i)*ones(size(x)), 'r');
    end
  elseif (nargin == 2)
    z_pos = axes;
    if (iscell(pts))
      for i=1:length(pts)
        plot3(pts{i}(:,1),pts{i}(:,2),ones(size(pts{i},1),1)*z_pos(i), 'k');
      end
    else
      indxs = unique(pts(:,4));
      for i=1:length(indxs)
        currents = (pts(:,4) == indxs(i));
        plot3(pts(currents,1),pts(currents,2),ones(sum(currents), 1)*z_pos(i), 'k');
      end
    end
  else
    if (iscell(pts))
      for i=1:length(pts)
        plot3(pts{i}(:,1),pts{i}(:,2),pts{i}(:,3), 'k');
      end
    else
      plot3(pts(:,1),pts(:,2),pts(:,3), 'k');
    end
  end

  hold off;

  return;
end
