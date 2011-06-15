function reconstruct_egg(mymovie, opts)

  [nframes, imgsize] = size_data(mymovie.dic);

  centers = mymovie.markers.centers;
  axes_length = mymovie.markers.axes_length;
  orientations = mymovie.markers.orientations;

  orientations = align_orientations(orientations);
  ratios = axes_length(1,:) ./ axes_length(2, :);

  target_ratio = median(ratios);
  new_axes = axes_length;
  new_axes(2, :) = sqrt(axes_length(1,:) .* axes_length(2,:) ./ target_ratio);
  new_axes(1, :) = new_axes(2, :) * target_ratio;

  mymovie = parse_metadata(mymovie, opts);

  [~, indx] = max(new_axes(1, :));
  equat_axes = new_axes(:, indx);

  z_pos = sqrt(equat_axes(2).^2 - new_axes(2, :).^2);

  figure;
  hold on;
  for i=1:nframes
    %plot3(x,y,z_pos(i)*ones(size(x)));
    %[x,y] = draw_ellipse(centers(:, i), axes_length(:, i), orientations(i));
    %plot3(x*opts.pixel_size,y*opts.pixel_size,mymovie.metadata.z_pos(i)*ones(size(x)), 'r');
    %[x,y] = draw_ellipse(centers(:, i), new_axes(:, i), orientations(i));
    %plot3(x*opts.pixel_size,y*opts.pixel_size,mymovie.metadata.z_pos(i)*ones(size(x)), 'b');
    plot3(mymovie.dic.eggshell(i).carth(:,1)*opts.pixel_size,mymovie.dic.eggshell(i).carth(:,2)*opts.pixel_size,mymovie.metadata.z_pos(i)*ones(size(mymovie.dic.eggshell(i).carth(:,1))), 'g');
  end
  axis('equal')

  figure;
  hold on;
  for i=1:nframes
    [x,y] = draw_ellipse([0;0], axes_length(:, i), 0);
    plot3(x*opts.pixel_size,y*opts.pixel_size,mymovie.metadata.z_pos(i)*ones(size(x)), 'r');
    [x,y] = draw_ellipse([0;0], new_axes(:, i), 0);
    plot3(x*opts.pixel_size,y*opts.pixel_size,mymovie.metadata.z_pos(i)*ones(size(x)), 'b');
  end
  axis('equal')

  return;
  %keyboard
  
  figure;
  hold on;
  for i=1:nframes
    pts = mymovie.markers.eggshell(i).carth;
    mid_pts = project2midplane(pts, centers(:,i), equat_axes, orientations(i), z_pos(i));
    %[x,y] = draw_ellipse([0;0], new_axes(:, i), 0);
    plot3(mid_pts(:, 1), mid_pts(:,2),i*ones(size(mid_pts(:,1))));
    plot3(pts(:, 1), pts(:,2),i*ones(size(mid_pts(:,1))), 'r');
  end

  figure;
  hold on;
  for i=1:nframes
    [x,y] = draw_ellipse(centers(:, i), new_axes(:, i), orientations(i));
    pts = mymovie.markers.eggshell(i).carth;
    plot3(pts(:, 1), pts(:,2),z_pos(i)*ones(size(pts(:,1))), 'g');
    plot3(x,y,z_pos(i)*ones(size(x)));
  end

  return;
end
