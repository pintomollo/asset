function [flows, dist] = flow_properties(mymovie, opts, data_only)

  if (nargin == 1)
    opts = [];
    data_only = false;
  elseif (nargin == 2)
    data_only = false;
  end

  is_char = (ischar(mymovie));
  thresh = [];
  if (is_char)
    fname = mymovie;
    files = dir(fname);

    all_data = [];
    dist = [];
    for i=1:length(files)
      load(files(i).name);
      [tmp_flows, tmp_dist] = flow_properties(mymovie, opts, true);

      if (i == 1)
        all_data = tmp_flows;
        dist = tmp_dist;
      else

        if (size(all_data, 1) < size(tmp_flows,1))
          all_data = [all_data tmp_flows(1:size(all_data,1),:,:)];
        else
          all_data = [all_data(1:size(tmp_flows,1),:, :) tmp_flows];
        end
      end
    end

    [h,w,n] = size(all_data);
    ndists = length(dist);

    all_avg = NaN([h, w, ndists]);
    all_x = all_data;
    ngoods = zeros(ndists,1);

    for i=1:ndists
      all_x(:,:,[2*i-1 2*i]) = i;
      all_avg(:,:,i) = nanmean(all_data(:,:,[2*i-1 2*i]),3);
      ngoods(i) = sum(sum(~isnan(all_avg(:,:,i))));
    end

    min_dist = [0 dist(1:end-1)];
  else
    if (exist([mymovie.experiment 'flows.mat'], 'file'))
      display('Loading previously computed data');
      load([mymovie.experiment 'flows.mat']);

      min_dist = [0 dist(1:end-1)];
      ndists = length(dist);
    else
      %dist = [0.125 0.25 0.5 1 1.25 1.5 1.75 2 3 5 10];
      dist = [0.5 1 1.5 2 2.5 3 4 5 10 15];
      min_dist = [0 dist(1:end-1)];
      ndists = length(dist);
      flows = cell(ndists, 1);

      opts.spot_tracking.projection_type = 'range';

      for i = 1:ndists
        opts.spot_tracking.projection_dist = dist(i);
        opts.spot_tracking.projection_dist_min = min_dist(i);
        mymovie = measure_flow(mymovie, opts);
        flows{i} = display_flow(mymovie, opts);
      end

      save([mymovie.experiment 'flows.mat'], 'flows', 'dist');
    end

    t = get_manual_timing(mymovie, opts);
    [h,w] = size(flows{1});
    half = (h/2);

    all_avg = NaN([half, t(1), ndists]);
    all_data = NaN([half, t(1), 2*ndists]);
    all_x = all_data;

    ngoods = zeros(ndists,1);

    for i=1:ndists
      all_data(:,:,2*i-1) = flows{i}(half+1:end, 1:t(1));
      all_data(:,:,2*i) = -flows{i}(half:-1:1, 1:t(1));
      all_x(:,:,[2*i-1 2*i]) = i;
      all_avg(:,:,i) = nanmean(all_data(:,:,[2*i-1 2*i]),3);
      ngoods(i) = sum(sum(~isnan(all_avg(:,:,i))));
    end

    full_data = all_data;

    avg = nanmean(all_data, 3);
    indx = find(all(isnan(avg), 1), 1, 'last') + 1;

    if (isempty(indx))
      indx = 1;
    end

    all_data = all_data(:, indx:end, :);
    all_x = all_x(:, indx:end, :);
    all_avg = all_avg(:, indx:end, :);

    avg = nanmean(all_data, 3);
    indx = find(any(isnan(avg), 2), 1, 'first') - 1;

    if (isempty(indx))
      indx = size(all_data, 1);
    end

    all_data = all_data(1:indx, :, :);
    all_x = all_x(1:indx, :, :);
    all_avg = all_avg(1:indx, :, :);
  end

  if (data_only)
    flows = all_data;
    return;
  end

  [height,width,n] = size(all_data);
  ngoods = ngoods / (height*width);

  indx = find(ngoods>0.8, 1, 'first');
  scaling = inpaint_nans(all_avg(:,:,indx));
  scaling = median_mex(scaling, 5);

  all_rescaled = bsxfun(@times, all_data, 1./scaling);
  goods = (abs(scaling) >= 0.01);
  all_goods = bsxfun(@and, ~isnan(all_rescaled), goods);

  H = 15;
  norm_range = (dist(indx) + min_dist(indx))/2;

  %inits = [0.1:0.1:10];
  inits = dist(2:end-1);
  best_h = NaN(size(inits));
  errs = best_h;
  for i=1:length(inits)
    [h0, r] = nlinfit(all_x(all_goods), all_rescaled(all_goods), @predicted_flow, inits(i));
    best_h(i) = h0;
    errs(i) = mean(abs(r));
  end

  [err, bindx] = min(errs);
  h = best_h(bindx);

  avg_vals = NaN(ndists,3);
  avg_stds = NaN(ndists,3);
  all_vals = [];
  for i=1:ndists
    vals = all_data(:,:,[2*i-1 2*i]);
    vals = vals(:);
    all_vals = [all_vals, vals];
    avg_vals(i,1) = nanmean(vals);
    avg_stds(i,1) = nanstd(vals);
    vals = all_data(1:round(height/2),:,[2*i-1 2*i]);
    vals = vals(:);
    avg_vals(i,2) = nanmean(vals);
    avg_stds(i,2) = nanstd(vals);
    vals = all_data(:,round(width/2):end,[2*i-1 2*i]);
    vals = vals(:);
    avg_vals(i,3) = nanmean(vals);
    avg_stds(i,3) = nanstd(vals);
  end

  pos = (dist + min_dist) / 2;
  rel_vals = bsxfun(@mrdivide, avg_vals, avg_vals(indx, :));
  rel_stds = bsxfun(@mrdivide, avg_stds, avg_vals(indx, :));
  pred_vals = predicted_flow(h, 1:ndists);

  subplot(111);
  hold off;
  boxplot(all_vals, 'positions', pos, 'labels', pos);
  hold on;
  errorbar(pos, rel_vals(:,1), rel_stds(:,1), 'r');
  errorbar(pos, rel_vals(:,2), rel_stds(:,2), 'g');
  errorbar(pos, rel_vals(:,3), rel_stds(:,3), 'b');
  plot(pos, pred_vals, 'k');
  ylim([-2 3]);
  title(num2str(h))

  if (is_char)
    print('-dpng', ['PNG/' fname 'avg.png']);

    for i=1:length(files)
      load(files(i).name);
      t = get_manual_timing(mymovie, opts);

      opts.spot_tracking.projection_args = Inf;
      mymovie = measure_flow(mymovie, opts);
      f0 = display_flow(mymovie, opts);
      
      opts.spot_tracking.projection_type = 'gaussian';
      opts.spot_tracking.projection_dist = h;
      opts.spot_tracking.projection_args = h;
      mymovie = measure_flow(mymovie, opts);
      f1 = display_flow(mymovie, opts);

      half = (size(f0,1)/2);

      s0 = (f0(half+1:end, 1:t(1)) -f0(half:-1:1, 1:t(1)))/2;
      s1 = (f1(half+1:end, 1:t(1)) -f1(half:-1:1, 1:t(1)))/2;

      subplot(121);
      hold off;
      imagesc(f0, [-0.2 0.2])
      subplot(122);
      hold off;
      imagesc(f1, [-0.2 0.2])
      print('-dpng', ['PNG/' mymovie.experiment 'flow-all.png']);

      subplot(121);
      hold off;
      imagesc(s0, [-0.2 0.2])
      subplot(122);
      hold off;
      imagesc(s1, [-0.2 0.2])
      print('-dpng', ['PNG/' mymovie.experiment 'half-all.png']);

    end
  else
    print('-dpng', ['PNG/' mymovie.experiment 'avg.png']);

    opts.spot_tracking.projection_args = Inf;
    mymovie = measure_flow(mymovie, opts);
    f0 = display_flow(mymovie, opts);
    
    opts.spot_tracking.projection_type = 'gaussian';
    opts.spot_tracking.projection_dist = h;
    opts.spot_tracking.projection_args = h;
    mymovie = measure_flow(mymovie, opts);
    f1 = display_flow(mymovie, opts);

    s0 = (f0(half+1:end, 1:t(1)) -f0(half:-1:1, 1:t(1)))/2;
    s1 = (f1(half+1:end, 1:t(1)) -f1(half:-1:1, 1:t(1)))/2;

    subplot(121);
    hold off;
    imagesc(f0, [-0.2 0.2])
    subplot(122);
    hold off;
    imagesc(f1, [-0.2 0.2])
    print('-dpng', ['PNG/' mymovie.experiment 'flow.png']);

    subplot(121);
    hold off;
    imagesc(s0, [-0.2 0.2])
    subplot(122);
    hold off;
    imagesc(s1, [-0.2 0.2])
    print('-dpng', ['PNG/' mymovie.experiment 'half.png']);
  end

  %inits = [0.1:0.1:10];
  %inits = dist;
  inits = dist(2:end-1);
  best_h = NaN(length(inits), 2);
  errs = NaN(length(inits), 1);
  for i=1:length(inits)
    [h0, r] = nlinfit(all_x(all_goods), all_rescaled(all_goods), @predicted_flow, [inits(i) 1]);
    best_h(i, :) = h0;
    errs(i) = mean(abs(r));
  end

  [err, bindx] = min(errs);
  h = best_h(bindx, :);

  avg_vals = NaN(ndists,3);
  avg_stds = NaN(ndists,3);
  all_vals = [];
  for i=1:ndists
    vals = all_data(:,:,[2*i-1 2*i]);
    vals = vals(:);
    all_vals = [all_vals, vals];
    avg_vals(i,1) = nanmean(vals);
    avg_stds(i,1) = nanstd(vals);
    vals = all_data(1:round(height/2),:,[2*i-1 2*i]);
    vals = vals(:);
    avg_vals(i,2) = nanmean(vals);
    avg_stds(i,2) = nanstd(vals);
    vals = all_data(:,round(width/2):end,[2*i-1 2*i]);
    vals = vals(:);
    avg_vals(i,3) = nanmean(vals);
    avg_stds(i,3) = nanstd(vals);
  end

  pos = (dist + min_dist) / 2;
  rel_vals = bsxfun(@mrdivide, avg_vals, avg_vals(indx, :));
  rel_stds = bsxfun(@mrdivide, avg_stds, avg_vals(indx, :));
  pred_vals = predicted_flow(h, 1:ndists);

  subplot(111);
  hold off;
  boxplot(all_vals, 'positions', pos, 'labels', pos);
  hold on;
  errorbar(pos, rel_vals(:,1), rel_stds(:,1), 'r');
  errorbar(pos, rel_vals(:,2), rel_stds(:,2), 'g');
  errorbar(pos, rel_vals(:,3), rel_stds(:,3), 'b');
  plot(pos, pred_vals, 'k');
  ylim([-2 3]);
  title(num2str(h))

  if (is_char)
    print('-dpng', ['PNG/' fname 'avg_2p.png']);

    for i=1:length(files)
      load(files(i).name);
      t = get_manual_timing(mymovie, opts);

      opts.spot_tracking.projection_args = Inf;
      mymovie = measure_flow(mymovie, opts);
      f0 = display_flow(mymovie, opts);
      
      opts.spot_tracking.projection_type = 'gaussian';
      opts.spot_tracking.projection_dist = h(1);
      opts.spot_tracking.projection_args = h(1);
      mymovie = measure_flow(mymovie, opts);
      f1 = display_flow(mymovie, opts);

      half = (size(f0,1)/2);

      s0 = (f0(half+1:end, 1:t(1)) -f0(half:-1:1, 1:t(1)))/2;
      s1 = (f1(half+1:end, 1:t(1)) -f1(half:-1:1, 1:t(1)))/2;

      subplot(121);
      hold off;
      imagesc(f0, [-0.2 0.2])
      subplot(122);
      hold off;
      imagesc(f1, [-0.2 0.2])
      print('-dpng', ['PNG/' mymovie.experiment 'flow-all_2p.png']);

      subplot(121);
      hold off;
      imagesc(s0, [-0.2 0.2])
      subplot(122);
      hold off;
      imagesc(s1, [-0.2 0.2])
      print('-dpng', ['PNG/' mymovie.experiment 'half-all_2p.png']);

    end
  else
    print('-dpng', ['PNG/' mymovie.experiment 'avg_2p.png']);

    opts.spot_tracking.projection_args = Inf;
    mymovie = measure_flow(mymovie, opts);
    f0 = display_flow(mymovie, opts);
    
    opts.spot_tracking.projection_type = 'gaussian';
    opts.spot_tracking.projection_dist = h(1);
    opts.spot_tracking.projection_args = h(1);
    mymovie = measure_flow(mymovie, opts);
    f1 = display_flow(mymovie, opts);

    s0 = (f0(half+1:end, 1:t(1)) -f0(half:-1:1, 1:t(1)))/2;
    s1 = (f1(half+1:end, 1:t(1)) -f1(half:-1:1, 1:t(1)))/2;

    subplot(121);
    hold off;
    imagesc(f0, [-0.2 0.2])
    subplot(122);
    hold off;
    imagesc(f1, [-0.2 0.2])
    print('-dpng', ['PNG/' mymovie.experiment 'flow_2p.png']);

    subplot(121);
    hold off;
    imagesc(s0, [-0.2 0.2])
    subplot(122);
    hold off;
    imagesc(s1, [-0.2 0.2])
    print('-dpng', ['PNG/' mymovie.experiment 'half_2p.png']);
  end


  return;

  function tmp_vals = predicted_flow(curr_h, x)

    if (numel(curr_h) > 1)
      coef = curr_h(2);
      curr_h = curr_h(1);
    else
      coef = 1;
    end

    v0 = curr_h/(curr_h - norm_range);
    vH = -coef*0.75*v0*curr_h/(H-curr_h);

    curr_vals = NaN(ndists, 1);
    for j=1:ndists
      a = min_dist(j);
      b = dist(j);

      if (curr_h >= b)
        curr_vals(j) = v0*(1-(b+a)/(2*curr_h));
      elseif (curr_h <= a)
        curr_vals(j) = vH*(1-(3*H^2 - 3*H*(b+a) + (b^2+b*a+a^2))/(3*(H-curr_h)^2));
      else
        curr_vals(j) = v0*(1-(curr_h+a)/(2*curr_h)) + vH*(1-(3*H^2 - 3*H*(b+curr_h) + (b^2+b*curr_h+curr_h^2))/(3*(H-curr_h)^2));
      end
    end

    tmp_vals = curr_vals(x);

    return;
  end
end
