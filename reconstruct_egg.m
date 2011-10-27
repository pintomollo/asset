function mymovie = reconstruct_egg(mymovie, opts)

  %%%%%%%%% ERRORS IN MOVIES 8 13 16 17 18


  dist_thresh = 50;

  [nframes, imgsize] = size_data(mymovie.dic);

  centers = mymovie.dic.centers;
  axes_length = mymovie.dic.axes_length;
  orientations = mymovie.dic.orientations;
  orientations = align_orientations(orientations);

  mymovie = parse_metadata(mymovie, opts);
  real_z = mymovie.metadata.z_position;

oks = all(axes_length > 0);
  real_thresh = median(axes_length(1,oks)) / dist_thresh;

r = robustfit(axes_length(1,oks), axes_length(2,oks));
v1 = bsxfun(@minus, axes_length, [0; r(1)]);

v2 = [1; r(2)];

c2 = dot(v1, repmat(v2, 1, size(v1, 2)));
c = dot(v2, v2);
b = c2/c;
r3 = robustfit([real_z(oks).^2; real_z(oks)].', b(oks).');
new_y = r3(2)*real_z.^2 + r3(3)*real_z+r3(1);

dist = abs(b - new_y);
goods = (abs(dist) < real_thresh);


  pts = cell(nframes, 1);
  all_coords = zeros(0,3);
  for i=1:nframes
    tmp_pts = mymovie.dic.eggshell(i).carth;
    if (~isempty(tmp_pts) & goods(i))
      npts = size(tmp_pts, 1);
      tmp_pts = tmp_pts * opts.pixel_size;
      all_coords = [all_coords; [tmp_pts ones(npts, 1)*real_z(i)]];
      pts{i} = tmp_pts;
    end
  end

  [fit_center, fit_axes, fit_angle] = ellipsoid_fit(all_coords);
  orient = atan2(-fit_angle(2,1), fit_angle(1,1));
  orient = orient + 2*pi*(orient < 0);

  z_coef = sqrt(1 - ((real_z - fit_center(3)).^2 / fit_axes(3).^2));
  new_axes = fit_axes(1:2, 1) * z_coef;

%  thresh = pi/16;
%  oks = (abs(orientations - orient) < thresh) & goods;

%  ratios = axes_length(1,:) ./ axes_length(2, :);
%  obj_ratio = fit_axes(1) / fit_axes(2);
%  thresh = 0.1;
%  oks = oks & (abs(ratios - obj_ratio) < thresh);

%  obj_center = fit_center(1:2);
%  dist = sqrt(sum(bsxfun(@minus, centers * opts.pixel_size, obj_center).^2, 1));
%  thresh = 20;
%  oks = oks & (dist < thresh);

%  dist = NaN(nframes, 1);
  npts = 128;
  dists = NaN(npts, nframes);
  pos = [0:npts-1] * (2*pi / npts);

%  figure;hold on;
  for i=1:nframes
    %pts = mymovie.dic.eggshell(i).carth;
    if (goods(i))
      tmp_pts = pts{i};
      npts = size(tmp_pts, 1);
      %pts = pts * opts.pixel_size;
      %plot3(pts(:,1), pts(:,2),ones(npts, 1)*real_z(i), 'b');

      [ell] = draw_ellipse(fit_center(1:2), new_axes(:, i), orient);
      nell = size(ell, 1);
      %plot3(ell(:,1), ell(:,2), ones(nell, 1)*real_z(i), 'r');

      ell_pts = carth2elliptic(tmp_pts, fit_center(1:2), new_axes(:, i), orient, 'radial');

      %dist(i) = mean(abs(ell_pts(:,2) - 1));
      [~, dists(:, i)] = interp_elliptic(ell_pts, pos);

    end
  end

  nclusters = 4;


  %%keyboard

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADDD NCLUSTER CHOOSING: E.G. JUMP METHOD

  clusters = zeros(1, nframes);
  correl = corr(dists(:, goods)); 
  clusters(goods) = kmeans(correl, nclusters, 'emptyaction' ,'drop', 'replicates', 5);
  %[~, indxs] = sort(clusters);
  
  cluster_count = zeros(nclusters, 1);
  for i=1:nclusters
    clusters_count(i) = sum(clusters == i);
  end 
  [~, best_cluster] = max(clusters_count);
  
  all_coords = zeros(0,3);
  for i=1:nframes
    if (clusters(i) == best_cluster)
      tmp_pts = pts{i};
      npts = size(tmp_pts, 1);
      %pts = pts * opts.pixel_size;
      all_coords = [all_coords; [tmp_pts ones(npts, 1)*real_z(i)]];
    end
  end

  [fit_center, fit_axes, fit_angle] = ellipsoid_fit(all_coords);
  orient = atan2(-fit_angle(2,1), fit_angle(1,1));
  orient = orient + 2*pi*(orient < 0);

  z_coef = sqrt(1 - ((real_z - fit_center(3)).^2 / fit_axes(3).^2));
  new_axes = fit_axes(1:2, 1) * z_coef;

  if (true)

  colors = ['mcyg'];
  colors(best_cluster) = 'r';

  corient = cos(-orient); 
  sorient = sin(-orient); 
  all_coords = zeros(0,3);

  figure;hold on;
  for i=1:nframes
    %pts = mymovie.dic.eggshell(i).carth;
    if (goods(i))
      tmp_pts = pts{i};
      npts = size(tmp_pts, 1);

      if (true)
      %pts = pts * opts.pixel_size;
      plot3(tmp_pts(:,1), tmp_pts(:,2),ones(npts, 1)*real_z(i), colors(clusters(i)));

      [ell] = draw_ellipse(fit_center(1:2), new_axes(:, i), orient);
      nell = size(ell, 1);
      plot3(ell(:,1), ell(:,2), ones(nell, 1)*real_z(i), 'k');
      end

      %ell_pts = carth2elliptic(pts, fit_center(1:2), new_axes(:, i), orient, 'radial');
      %dist(i) = mean(abs(ell_pts(:,2) - 1));
      if (clusters(i) == best_cluster)
        tmp_pts = [tmp_pts ones(npts, 1)*real_z(i)];
        tmp_pts = bsxfun(@minus, tmp_pts, fit_center.');

        x = (tmp_pts(:,1)*corient + tmp_pts(:,2)*sorient);
        y = (tmp_pts(:,1)*sorient - tmp_pts(:,2)*corient);

        all_coords = [all_coords; [x y tmp_pts(:,3)]];
      end
    end
  end

  %figure;plot3(all_coords(:, 1), all_coords(:, 2), all_coords(:, 3))

  %keyboard

  [fit_center2, fit_axes2, fit_angle2] = ellipsoid_fit(all_coords, 1);
  orient2 = atan2(-fit_angle2(2,1), fit_angle2(1,1));
  orient2 = orient2 + 2*pi*(orient2 < 0);

  [fit_axes fit_axes2]
  fit_axes = fit_axes2;

  %[fit_center fit_center2]
  %[orient orient2]

  end

  if (fit_axes(2) > fit_axes(1))
    fit_axes = fit_axes([2 1 3]);
  end

  mymovie.metadata.axes_length = fit_axes;

  if (false)
  z_coef = sqrt(1 - ((real_z - fit_center(3)).^2 / fit_axes(3).^2));
  new_axes = fit_axes(1:2, 1) * z_coef;

  i = round(nframes / 2);
%for i=1:nframes

  %if (~goods(i))
  %  continue;
  %end

  nimg = i;
  imshow(imnorm(double(load_data(mymovie.dic, nimg))));
  hold on;
  draw_ellipse(mymovie.dic.centers(:, nimg), mymovie.dic.axes_length(:, nimg), mymovie.dic.orientations(1, nimg));
  draw_ellipse(fit_center(1:2) / opts.pixel_size, new_axes(:, nimg) / opts.pixel_size, orient, 'r');

  title(real_z(nimg) - fit_center(3));
  drawnow
  saveas(gca, [mymovie.experiment '-Zpos-' num2str(nimg) '.jpg']);
  hold off;

  end
%end

%  z_pos = unique(all_coords(:,3)).';
%  figure;scatter3(all_coords(:,1), all_coords(:,2), all_coords(:,3));
%  hold on;

%  z_coef = sqrt(1 - ((z_pos.^2) / (fit_axes2(3).^2)));
%  new_axes = fit_axes2(1:2, 1) * z_coef;

%  for i=length(z_pos)
%
%      [ell] = draw_ellipse(fit_center2(1:2), new_axes(:, i), orient2);
%      nell = size(ell, 1);
%      plot3(ell(:,1), ell(:,2), ones(nell, 1)*z_pos(i), 'k');
%
%  end

  %figure;
  %subplot(2,2,1)
  %plot(hypot(centers(1,:) - fit_center(1), centers(2,:) - fit_center(2)))
  %subplot(2,2,2)
  %plot(abs((axes_length(1,:) ./ axes_length(2, :)) - (fit_axes(1) / fit_axes(1))));
  %subplot(2,2,3)
  %plot(abs(orientations - orient));
  %subplot(2,2,4)
  %plot(dist);


%  keyboard

  return;

  goods = true(nframes, 1);
  paths = cell(nframes, 1);


  %figure;scatter3(real_z, axes_length(1,:), axes_length(2,:))
  %figure;scatter(axes_length(1,:), axes_length(2,:))

oks = all(axes_length > 0);
  real_thresh = median(axes_length(1,oks)) / dist_thresh;

r = robustfit(axes_length(1,oks), axes_length(2,oks));
v1 = bsxfun(@minus, axes_length, [0; r(1)]);

v2 = [1; r(2)];

c2 = dot(v1, repmat(v2, 1, size(v1, 2)));
c = dot(v2, v2);
b = c2/c;
r3 = robustfit([real_z(oks).^2; real_z(oks)].', b(oks).');
new_y = r3(2)*real_z.^2 + r3(3)*real_z+r3(1);

dist = abs(b - new_y);
goods = (abs(dist) < real_thresh);

z0 = -r3(3)/(2*r3(2));
real_axes = axes_length * opts.pixel_size;


  obj_orient = median(orientations(goods));
  thresh = pi/16;
  oks = (abs(orientations - obj_orient) < thresh) & goods;

  ratios = axes_length(1,:) ./ axes_length(2, :);
  obj_ratio = median(ratios(oks));
  thresh = 0.1;
  oks = oks & (abs(ratios - obj_ratio) < thresh);

  obj_center = mean(centers(:, oks), 2);
  dist = sqrt(sum(bsxfun(@minus, centers, obj_center).^2, 1));
  thresh = 20;
  oks = oks & (dist < thresh);



goods = oks;


ratios = real_axes(1,:) ./ real_axes(2, :);
obj_ratio = median(ratios(goods));

i = 1;
coefs = zeros(3, 1);
prev_coefs = coefs + 1;
while(any(coefs ~= prev_coefs))
  i = i + 1;
  prev_coefs = coefs;
  [coefs, cntrs] = estimate_axes(real_axes(:, goods), real_z(goods) - z0, obj_ratio);
  obj_ratio = coefs(1) / coefs(2);
  z0 = z0 + cntrs(3);

  if (i > 100)
    break;
  end
end

mymovie.metadata.axes_length = coefs;
mymovie.metadata.z_center = z0;

frames = randperm(sum(goods));
good_frames = [1:nframes];
good_frames = good_frames(goods);
frames = good_frames(frames);
rel_z = real_z - z0;

z_coef = sqrt(1 - (rel_z.^2 / (coefs(3)^2)));
new_axes = coefs(1:2, 1) * z_coef;
new_axes = new_axes / opts.pixel_size;

for i=frames
  nimg = i;
  imshow(imnorm(double(load_data(mymovie.dic, nimg))));
  hold on;
  draw_ellipse(mymovie.dic.centers(:, nimg), mymovie.dic.axes_length(:, nimg), mymovie.dic.orientations(1, nimg));
  draw_ellipse(mymovie.dic.centers(:, nimg), new_axes(:, nimg), mymovie.dic.orientations(1, nimg), 'r');

  title(rel_z(nimg));
  drawnow
  saveas(gca, [mymovie.experiment '-Zpos-' num2str(nimg) '.jpg']);
  hold off;
end

return;

%[c1, a1, o1] = fit_ellipse([real_axes(1,goods) -real_axes(1,goods)], [real_z(goods) real_z(goods)] - z0);
%[c2, a2, o2] = fit_ellipse([real_axes(2,goods) -real_axes(2,goods)], [real_z(goods) real_z(goods)] - z0);
%[c3, a3, o3] = fit_ellipse([[real_axes(1,goods) -real_axes(1,goods)], [real_axes(2,goods) -real_axes(2,goods)]*obj_ratio], [real_z(goods) real_z(goods) real_z(goods) real_z(goods)]);

%z_coef = sqrt(1 - ((real_z - c3(2)).^2 / (a3(2)^2)));
%z_proj = bsxfun(@rdivide, real_axes, z_coef);

%nelems = sum(goods);
%[c1, a1, o1] = fit_ellipse([z_proj(1, goods) zeros(1, nelems)], [zeros(1, nelems) z_proj(2, goods)]);

keyboard

%rel_z = real_z - 

figure;
plot(real_z, b);
hold on;
plot(real_z, new_y, 'g');
plot(real_z(goods), b(goods), 'r');

%  keyboard

%  p1 = polyfit(real_z, axes_length(1,:), 2);
%  p2 = polyfit(real_z, axes_length(2,:), 2);
%  y1 = polyval(p1, real_z);
%  y2 = polyval(p2, real_z);

%  figure;
%  subplot(2,1,1)
%  plot(real_z, axes_length(1, :));
%  hold on;
%  plot(real_z, y1, 'r');
%  subplot(2,1,2)
%  plot(real_z, axes_length(2, :));
%  hold on;
%  plot(real_z, y2, 'r');


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

  oks = goods;

  obj_orient = median(orientations(goods));
  thresh = pi/16;
  oks = (abs(orientations - obj_orient) < thresh);

  ratios = axes_length(1,:) ./ axes_length(2, :);
  obj_ratio = median(ratios(oks));
  thresh = 0.1;
  oks = oks & (abs(ratios - obj_ratio) < thresh);

  obj_center = mean(centers(:, oks), 2);
  dist = sqrt(sum(bsxfun(@minus, centers, obj_center).^2, 1));
  thresh = 20;
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
  relative_z = real_z - c(2);
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
  subplot(2,1,1)
  plot(axes_length(1, :), relative_z);
  subplot(2,1,2)
  plot(axes_length(2, :), relative_z);
  figure;
  plot3(centers(1,:), centers(2,:), relative_z);
  hold on;
  plot3(centers(1,oks), centers(2,oks), relative_z(oks), 'r');
  figure;
  plot(ratios, relative_z);
  hold on;
  plot(ratios(oks), relative_z(oks), 'r');
  figure;
  plot(orientations, relative_z);
  hold on;
  plot(orientations(oks), relative_z(oks), 'r');

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

    rescale = sqrt(1 - (relative_z(i)/coefs(3))^2);
    %[centers(:,i), axes_length(:,i), orientations(1,i)] = fit_ellipse(pts);

    %plot3(x,y,z_pos(i)*ones(size(x)));
    [x,y] = draw_ellipse([0;0], coefs(1:2)*rescale, 0);
    plot3(x,y, relative_z(i)*ones(size(x)), 'r');

    %[x,y] = draw_ellipse(centers(:, i), new_axes(:, i), orientations(i));
    %plot3(x*opts.pixel_size,y*opts.pixel_size,mymovie.metadata.z_pos(i)*ones(size(x)), 'b');
    plot3(pts(:,1), pts(:,2), relative_z(i)*ones(size(pts(:,1))), 'g');
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


function [axes_size, centers] = estimate_axes(axes_length, z_pos, axes_ratio)

[c3, a3, o3] = fit_ellipse([[axes_length(1,:), -axes_length(1,:)], [axes_length(2,:), -axes_length(2,:)]*axes_ratio], repmat(z_pos, 1, 4));

z_coef = sqrt(1 - ((z_pos - c3(2)).^2 / (a3(2)^2)));
z_proj = bsxfun(@rdivide, axes_length, z_coef);

nelems = length(z_pos);
[c1, a1, o1] = fit_ellipse([z_proj(1, :) zeros(1, nelems) -z_proj(1,:) zeros(1,nelems)], [zeros(1, nelems) z_proj(2, :) zeros(1,nelems) -z_proj(2,:)]);

  axes_size = [a1; a3(2)];
  centers = [c1; c3(2)];

  return;
end
