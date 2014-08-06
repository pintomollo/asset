function mymovie = piv_flows(mymovie, opts)

  obj_dir = fullfile('PIV', mymovie.experiment);
  if (~exist('PIV') | ~exist(obj_dir))
    warning('No relevant data found')

    return;
  end

  files = dir(fullfile(obj_dir, 'seq_*_PIV3_disp.txt'));
  nframes = length(mymovie.data.cortex);
  ncols = 16;

  proj_dist = opts.spot_tracking.projection_dist / opts.pixel_size;
  pos_bin = opts.spot_tracking.projection_bin_size;
  frame_bin = opts.spot_tracking.projection_frames;
  if (isfield(opts.spot_tracking, 'projection_dist_min'))
    proj_min = opts.spot_tracking.projection_dist_min / opts.pixel_size;
  else
    proj_min = 0;
  end
  if (isfield(opts.spot_tracking, 'projection_args'))
    proj_args = opts.spot_tracking.projection_args;
  else
    proj_args = [];
  end
  if (strncmp(opts.spot_tracking.projection_type, 'gaussian', 8))
    if (isempty(proj_args))
      nullcline = 1.4;
    else
      nullcline = proj_args(1);
      proj_args = proj_args(2:end);
    end
    g_factor = 0.85;
  end

  if length(proj_args) > 1
    resample = proj_args(1);
    remove_outside = proj_args(2);
  else
    resample = 4;
    remove_outside = true;
  end

  pivs = cell(nframes, 1);

  for i=1:length(files)
    fname = fullfile(obj_dir, files(i).name);
    data = textread(fname, '%n');
    frame = regexp(files(i).name, 'seq_(\d+)_PIV', 'tokens');

    if (~isempty(frame))
      frame = str2double(frame{1})+1;
      data = reshape(data, ncols, []).';
      npts = sqrt(size(data, 1));
      data = reshape(data, [npts npts, ncols]);
      pivs{frame} = imresize(data, resample);
    end
  end

  edges = [0:pos_bin:65];
  edges = [-edges(end:-1:2) edges].';
  centers = edges(1:end-1) + pos_bin/2;
  edges(1) = -Inf;
  edges(end) = Inf;

  tmp_speed = NaN(size(centers));
  tmp_std = tmp_speed;

  all_speed = [];
  all_pos = [];
  all_frame = [];

  nucleus_center = (isfield(mymovie.data, 'nuclei') & ~empty_struct(mymovie.data.nuclei, 'carth'));
  if (nucleus_center)
    for i=1:nframes
      if (~isempty(mymovie.data.nuclei(i).carth))
        prev_center = [mymovie.data.nuclei(i).carth mymovie.data.nuclei(i).properties i];
        break;
      end
    end
  end

  nucleus_width = 1.25;

  all_dists = NaN(nframes, 2);

  init_shift = 1;

  for i=(init_shift+1):nframes
    nimg = i;

    flows = reshape(pivs{nimg}, [], ncols);

    cortex = mymovie.data.cortex(nimg).carth;
    ncortex = size(cortex, 1);

    if (~isempty(flows) & ncortex > 0)

      egg_angle = mymovie.data.orientations(nimg);
      if (nucleus_center)
        if (~isempty(mymovie.data.nuclei(nimg).carth))
          nucleus_center = [mymovie.data.nuclei(nimg).carth mymovie.data.nuclei(nimg).properties nimg];
        else
          nucleus_center = prev_center;
        end

        ells = carth2elliptic([cortex; nucleus_center(1:2)], mymovie.data.centers(:, nimg), [1;1], 0, 'radial');
        ell_nucleus = ells(end,:);
        ells = ells(1:end-1, :);

        dangle = asin(nucleus_width*nucleus_center(3)/ell_nucleus(2));
  
        egg_angle = align_orientations(egg_angle, ell_nucleus(1));

        good_ells = (ells(:,1) >= ell_nucleus(1) - dangle & ells(:,1) <= ell_nucleus(1) + dangle);

        if (ell_nucleus(1) - dangle < 0)
          good_ells = good_ells | (ells(:,1)-2*pi >= ell_nucleus(1) - dangle & ells(:,1)-2*pi <= ell_nucleus(1) + dangle);
        end
        if (ell_nucleus(1) + dangle > 2*pi)
          good_ells = good_ells | (ells(:,1)+2*pi >= ell_nucleus(1) - dangle & ells(:,1)+2*pi <= ell_nucleus(1) + dangle);
        end

        dist = sqrt(mymean(sum(bsxfun(@minus, cortex(good_ells,:), nucleus_center(1:2)).^2, 2), 1));
        all_dists(i,1) = dist;
        dist = dist / nucleus_center(3);
        %rads & some function ....
        all_dists(i,2) = dist;

        %dangle = dangle*dist;
        % Sigmoid between 2 and 4 roughly, 3./(1+exp(-5*(x-2.5)))+1

        dangle = dangle*(4./(1+exp(-5*(dist-2.5)))+1);

        if (dist <= 2)
          egg_angle = ell_nucleus(1);
        elseif (dist <= 3)
          dist = dist - 2;
          egg_angle = ell_nucleus(1)*(1-dist) + egg_angle*dist;
        end

        good_ells = (ells(:,1) >= ell_nucleus(1) - dangle & ells(:,1) <= ell_nucleus(1) + dangle);

        if (ell_nucleus(1) - dangle < 0)
          good_ells = good_ells | (ells(:,1)-2*pi >= ell_nucleus(1) - dangle & ells(:,1)-2*pi <= ell_nucleus(1) + dangle);
        end
        if (ell_nucleus(1) + dangle > 2*pi)
          good_ells = good_ells | (ells(:,1)+2*pi >= ell_nucleus(1) - dangle & ells(:,1)+2*pi <= ell_nucleus(1) + dangle);
        end


        prev_center = nucleus_center;
      end

      align_cortex = realign(cortex, [0;0], mymovie.data.centers(:, nimg), egg_angle);
      [pos, post_indxs, tot_dist] = carth2linear(align_cortex, true, opts);

      cortex = cortex(post_indxs, :);
      [perp] = perpendicular_sampling(cortex, opts);
      pos = pos * tot_dist;

      if (remove_outside)
        ins = inpolygon(flows(:,1), flows(:,2),cortex(:,1), cortex(:,2));
        flows = flows(ins,:);
      end

      speed = bsxfun(@times, flows(:,3:4), flows(:,5));
      pts = flows(:,1:2);
      dist = sqrt(bsxfun(@minus, cortex(:,1), pts(:,1).').^2 + ...
                  bsxfun(@minus, cortex(:,2), pts(:,2).').^2);
      dist(dist>3*proj_dist) = Inf;

      dist_thresh = nullcline / opts.pixel_size;

      weights = exp(-(dist.^2)/(2*((proj_dist/2)^2)));

      if (isfinite(dist_thresh))
        rescale = min(dist, [], 1);
        goods = (rescale < g_factor*dist_thresh);
        rescale(~goods) = 0;
        weights(:, ~goods) = 0;
        weights = bsxfun(@rdivide, weights, sum(weights, 2));
        rescale = bsxfun(@times, weights, rescale);
        rescale = sum(rescale, 2);
        goods = (rescale < g_factor*dist_thresh);
        scaling = dist_thresh ./ (dist_thresh - rescale);
        scaling(~goods) = NaN;
      else
        scaling = ones(size(dist, 1), 1);
      end

      weights = bsxfun(@rdivide, weights, sum(weights, 2));
      speed_x = sum(bsxfun(@times, weights, speed(:,1).'), 2);
      speed_y = sum(bsxfun(@times, weights, speed(:,2).'), 2);

      movement = [speed_x, speed_y];

      if (nucleus_center)
   %     keyboard
        movement(good_ells(post_indxs), :) = NaN;
      end

      speed = dot([perp(:, 2), -perp(:, 1)], movement, 2) .* scaling;

      all_speed = [all_speed; speed*opts.pixel_size];
      all_pos = [all_pos; pos*opts.pixel_size];
      all_frame = [all_frame; nimg*ones(size(speed))];

      currents = (all_frame > nimg-frame_bin);
      all_speed = all_speed(currents);
      all_pos = all_pos(currents);
      all_frame = all_frame(currents);

      curr_speed = tmp_speed;
      curr_std = tmp_std;
      [counts, groups] = histc(all_pos, edges);
      [curr_speed(counts~=0), curr_std(counts~=0), groups] = mymean(all_speed, 1, groups);

      flow(nimg).speed = [curr_speed, curr_std];
      flow(nimg).position = centers;
    else
      flow(nimg).speed = [];
      flow(nimg).position = [];
    end
  end

  mymovie.data.piv = flow;

  return;
end
