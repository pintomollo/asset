function mymovie = reconstruct_egg(mymovie, opts)

  [nframes, imgsize] = size_data(mymovie.dic);

  centers = mymovie.dic.centers;
  axes_length = mymovie.dic.axes_length;
  orientations = mymovie.dic.orientations;
  orientations = align_orientations(orientations);

  mymovie = parse_metadata(mymovie, opts);
  real_z = mymovie.metadata.z_position;

  goods = true(nframes, 1);
  paths = cell(nframes, 1);

  %figure;
  %hold on;
  for i=1:nframes
    pts = mymovie.dic.eggshell(i).carth;
    if (isempty(pts))
      goods(i) = false;
%      continue;
    end
    paths{i} = pts;
    %[centers(:,i), axes_length(:,i), orientations(1,i)] = fit_ellipse(pts);

    %plot3(x,y,z_pos(i)*ones(size(x)));
%    [x,y] = draw_ellipse(centers(:, i), axes_length(:, i), orientations(i));
%    plot3(x*opts.pixel_size,y*opts.pixel_size,real_z(i)*ones(size(x)), 'r');
    %[x,y] = draw_ellipse(centers(:, i), new_axes(:, i), orientations(i));
    %plot3(x*opts.pixel_size,y*opts.pixel_size,mymovie.metadata.z_pos(i)*ones(size(x)), 'b');
%    plot3(pts(:,1)*opts.pixel_size,pts(:,2)*opts.pixel_size, real_z(i)*ones(size(pts(:,1))), 'g');
  end

  centers(:, ~goods) = NaN;
  axes_length(:, ~goods) = NaN;
  orientations(1, ~goods) = NaN;

  obj_orient = median(orientations(goods));
  thresh = pi/32;
  oks = (abs(orientations - obj_orient) < thresh);

  ratios = axes_length(1,:) ./ axes_length(2, :);
  obj_ratio = median(ratios(oks));
  thresh = 0.05;
  oks = oks & (abs(ratios - obj_ratio) < thresh);

  obj_center = mean(centers(:, oks), 2);
  dist = sqrt(sum(bsxfun(@minus, centers, obj_center).^2, 1));
  thresh = 10;
  oks = oks & (dist < thresh);

  obj_center = mean(centers(:, oks), 2);
  obj_ratio = mean(ratios(oks));
  obj_orient = mean(orientations(goods));

  corient = cos(obj_orient); 
  sorient = sin(obj_orient); 

  all_pts = zeros(0,3);
  all_coords = zeros(0,3);

  for i=1:nframes
    if (~oks(i))
      continue;
    end
    pts = paths{i};
    pts = bsxfun(@minus, pts, obj_center.') * opts.pixel_size;

    pts(:,2) = -pts(:,2);

    x = (pts(:,1)*corient + pts(:,2)*sorient);
    y = (pts(:,1)*sorient - pts(:,2)*corient);
    paths{i} = [x -y];

    angl = atan2(y*obj_ratio, x);
    proj = sign(angl) .* sqrt(x.^2 + (obj_ratio^2) * (y.^2));
    z = real_z(i)*ones(size(proj));

    all_pts = [all_pts; [proj, proj/obj_ratio, z]];
    all_coords = [all_coords; [x y z]];
  end

  %figure;scatter(all_pts(:,1), all_pts(:,3));
  [c,a,o] = fit_ellipse(all_pts(:,[1 3]));
  all_coords(:,[1 3]) = bsxfun(@minus, all_coords(:,[1 3]), c.');
  real_z = real_z - c(2);
  %hold on;
  %draw_ellipse(c,a,o, 'r')

  %figure;scatter(all_pts(:,2), all_pts(:,3));
  [c,a,o] = fit_ellipse(all_pts(:,[2 3]));
  all_coords(:,2) = all_coords(:,2) - c(1);

  b = ones(size(all_coords(:,1)));
  coefs = all_coords.^2 \ b;
  coefs = 1 ./ sqrt(coefs);
  mymovie.metadata.axes_length = coefs;

  %hold on;
  %draw_ellipse(c,a,o, 'r')

  if (true)
  figure;
  plot3(centers(1,:), centers(2,:), real_z);
  hold on;
  plot3(centers(1,oks), centers(2,oks), real_z(oks), 'r');
  figure;
  plot(ratios, real_z);
  hold on;
  plot(ratios(oks), real_z(oks), 'r');
  figure;
  plot(orientations, real_z);
  hold on;
  plot(orientations(oks), real_z(oks), 'r');

  figure;
  hold on;
  for i=1:nframes

    %pts = mymovie.dic.eggshell(i).carth;
    if (~oks(i))
      %pts = paths{i};
      %if (~isempty(pts))
      %  plot3(pts(:,1), pts(:,2), real_z(i)*ones(size(pts(:,1))), 'y');
      %end
      continue;
    end
    pts = paths{i};

    rescale = sqrt(1 - (real_z(i)/coefs(3))^2);
    %[centers(:,i), axes_length(:,i), orientations(1,i)] = fit_ellipse(pts);

    %plot3(x,y,z_pos(i)*ones(size(x)));
    [x,y] = draw_ellipse([0;0], coefs(1:2)*rescale, 0);
    plot3(x,y, real_z(i)*ones(size(x)), 'r');

    %[x,y] = draw_ellipse(centers(:, i), new_axes(:, i), orientations(i));
    %plot3(x*opts.pixel_size,y*opts.pixel_size,mymovie.metadata.z_pos(i)*ones(size(x)), 'b');
    plot3(pts(:,1), pts(:,2), real_z(i)*ones(size(pts(:,1))), 'g');
  end
  end

  return;

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
  real_z = mymovie.metadata.z_pos;
  npts = length(real_z);

  lower = [diff(real_z)~=0 0];

  x_pos = [0:npts-1];
  derivatives = (39*(real_z(:,5:end-2) - real_z(:,3:end-4)) + 12*(real_z(:,6:end-1) - real_z(:,2:end-5)) -5*(real_z(:, 7:end) - real_z(:, 1:end-6))) / 96;
  y_valids = (derivatives > 0 & ~isnan(derivatives));
  x_valids = logical([0 0 0 y_valids 0 0 0]);

  y_deriv = log(derivatives(y_valids));
  x_deriv = x_pos(x_valids).';

  expfunc = @(p,x)(p(1)*(1-exp(-p(2)*x)) + p(3));

  params = [x_deriv ones(size(x_deriv))] \ y_deriv(:);
  params = [-exp(params(2)) / params(1), -params(1), real_z(1)];

  valids = ~isnan(real_z) & lower;
  better_params = lsqcurvefit(expfunc,params,x_pos(valids), real_z(valids));
  relative_z = real_z - expfunc(better_params, x_pos);

  mymovie.metadata.z_rel = relative_z;

  figure;scatter(prod(new_axes), relative_z);
  figure;scatter(axes_length(1,:), relative_z);

  figure;plot(x_pos, real_z);
  hold on;plot(x_pos, expfunc(params, x_pos), 'r');
  plot(x_pos, expfunc(better_params, x_pos), 'g');

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
    h1 = subplot(211);
  hold on;
    [x,y] = draw_ellipse([0;0], axes_length(:, i), 0);
    plot3(x*opts.pixel_size,y*opts.pixel_size,relative_z(i)*ones(size(x)), 'r');
    h2 = subplot(212);
  hold on;
    [x,y] = draw_ellipse([0;0], new_axes(:, i), 0);
    plot3(x*opts.pixel_size,y*opts.pixel_size,relative_z(i)*ones(size(x)), 'b');
  end
  axis(h1, 'equal')
  axis(h2, 'equal')
  %linkaxes([h1 h2])

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
