function flow_nate(niter, use_timing, pos)

% Optimized alignement of flow measurements:
% p = [260; 337; 77; 337; 1; 337]



  aligned_times = [];
  if (nargin == 0)
    niter = 1;
    use_timing = false;
  elseif (nargin == 1) 
    if (islogical(niter))
      use_timing = niter;
      niter = 1;
    else
      use_timing = false;
    end
  elseif (isnumeric(use_timing) && numel(use_timing)>1)
    aligned_times = use_timing;
    use_timing = true;
  end

  if (iscell(niter))
    signals = niter;
    signals_full = stack_images(signals, aligned_times);
    nfiles = length(signals);
  else

    %files = dir('749-*_.mat');
    files = [dir('749-0*_.mat'); dir('749-2*_.mat')];
    nfiles = length(files);
    res = cell(nfiles, 3);
    signals = cell(nfiles, 1);
    times = NaN(nfiles, 1);
    signals_full = cell(niter, 3);

    for i=1:nfiles
      load(files(i).name);
      display(mymovie.experiment);

      %opts.spot_tracking.projection_bin_size = 2;
      
      opts.spot_tracking.projection_args = 2.49;
      opts.spot_tracking.projection_dist = 2.49/3;
      mymovie = measure_flow(mymovie, opts);
      time = get_manual_timing(mymovie, opts);
      [tmp, pos] = display_flow(mymovie, opts);
      tmp = tmp.';

      times(i) = time(2);

      %imagesc(tmp);
      %set(gca, 'XTickLabel', pos(get(gca, 'XTick')));
      %colorbar;
      %print('-dpng', ['PNG/' mymovie.experiment 'flow.png']);

      signals{i} = tmp;
    end

    if (use_timing && ~isempty(aligned_times))
      times = aligned_times;
    end

    prev_best = Inf;
    for j=1:3
      for i=1:niter
        if (use_timing)
          [signals_full{i,1}, signals_full{i,2}, signals_full{i,3}] = simultaneous_registration(signals, times - min(times) + 1);
        else
          [signals_full{i,1}, signals_full{i,2}, signals_full{i,3}] = simultaneous_registration(signals);
        end

      end

      uuid = num2str(now + cputime);
      display(num2str(uuid))
      save(['aligned_flow_' num2str(use_timing) '_' uuid '.mat'], 'signals_full', 'signals');

      vals = cat(1, signals_full{:,3});
      [v, indx] = min(vals);

      if (v < prev_best)
        best_signals = signals_full{indx, 1};
        prev_best = v;
      end

    %return;
    end
    signals_full = best_signals;
  end

%  signals2 = simultaneous_registration(signals2);
  signals = (signals_full(:,(end/2)+1:end,:) - signals_full(:,end/2:-1:1,:)) / 2;

%  imagesc(mymean(signals,3));
%  set(gca, 'XTickLabel', pos(get(gca, 'XTick')));
%  colorbar;
%  print('-dpng', ['PNG/flow_avg_half.png']);

  imagesc(mymean(signals,3));
  set(gca, 'XTickLabel', pos(get(gca, 'XTick')));
  colorbar;
  print('-dpng', ['PNG/flow_avg.png']);

  %{
  valids = (sum(~isnan(signals), 3) <= 4);
  signals(repmat(valids, [1 1 2*nfiles])) = NaN;

  signals = padarray(signals, [10 3 0], 0);
  p = pos((end/2)+1:end).';
  t = [0:size(signals, 1)-1].';
  p = [[-2:0], p, [67:69]];

  [X,Y] = meshgrid(p, t);
  X = repmat(X, [1 1 2*nfiles]);
  Y = repmat(Y, [1 1 2*nfiles]);

  goods = any(~valids, 2);

  inter = gridfit(X(:), Y(:), signals(:), p, t, 'smooth', [25 500]);
  inter = interp2(X(:,:,1), Y(:,:,1), inter, [1:66], t(goods));
  inter = padarray(inter, [1 1], 0);
  %}

  valids = (sum(~isnan(signals), 3) <= 2);
  signals(repmat(valids, [1 1 nfiles])) = 0;

  signals = padarray(signals, [10 3 0], 0);
  p = pos((end/2)+1:end);
  t = [0:size(signals, 1)-1].';
  p = [[-2:0].'; p; [67:69].'];

  [X,Y] = meshgrid(p, t);
  X = repmat(X, [1 1 nfiles]);
  Y = repmat(Y, [1 1 nfiles]);

  goods = any(~valids, 2);
  range_t = t(goods);
  range_t = range_t([1 end]);

  inter2 = gridfit(X(:), Y(:), signals(:), p, t, 'smooth', [25 500]);
  flow_full = interp2(X(:,:,1), Y(:,:,1), inter2, [1:66], t(goods));
  flow_full = padarray(flow_full, [10 1], 0);

  flow_small = interp2(X(:,:,1), Y(:,:,1), inter2, [1:66], [range_t(1):15:range_t(end)].');
  flow_small = padarray(flow_small, [1 1], 0);

  %imagesc(inter);
  %set(gca, 'XTickLabel', pos(get(gca, 'XTick')));
  %colorbar;
  %print('-dpng', ['PNG/flow_avg_half_smooth.png']);

  imagesc(flow_full);
  set(gca, 'XTickLabel', pos(get(gca, 'XTick')));
  colorbar;
  print('-dpng', ['PNG/flow_avg_smooth.png']);

  signals = signals_full(:,(end/2)+1:end,:);

  valids = (sum(~isnan(signals), 3) <= 2);
  signals(repmat(valids, [1 1 nfiles])) = 0;

  signals = padarray(signals, [10 3 0], 0);
  p = pos((end/2)+1:end);
  t = [0:size(signals, 1)-1].';
  p = [[-2:0].'; p; [67:69].'];

  [X,Y] = meshgrid(p, t);
  X = repmat(X, [1 1 nfiles]);
  Y = repmat(Y, [1 1 nfiles]);

  goods = any(~valids, 2);
  range_t = t(goods);
  range_t = range_t([1 end]);

  inter2 = gridfit(X(:), Y(:), signals(:), p, t, 'smooth', [25 500]);
  flow_full = interp2(X(:,:,1), Y(:,:,1), inter2, [1:66], t(goods));
  flow_full = padarray(flow_full, [10 1], 0);

  flow_small = interp2(X(:,:,1), Y(:,:,1), inter2, [1:66], [range_t(1):15:range_t(end)].');
  flow_small = padarray(flow_small, [1 1], 0);

  imagesc(flow_full);
  set(gca, 'XTickLabel', pos(get(gca, 'XTick')));
  colorbar;
  print('-dpng', ['PNG/flow_right_smooth.png']);

  signals = - signals_full(:,end/2:-1:1,:);

  valids = (sum(~isnan(signals), 3) <= 2);
  signals(repmat(valids, [1 1 nfiles])) = 0;

  signals = padarray(signals, [10 3 0], 0);
  p = pos((end/2)+1:end);
  t = [0:size(signals, 1)-1].';
  p = [[-2:0].'; p; [67:69].'];

  [X,Y] = meshgrid(p, t);
  X = repmat(X, [1 1 nfiles]);
  Y = repmat(Y, [1 1 nfiles]);

  goods = any(~valids, 2);
  range_t = t(goods);
  range_t = range_t([1 end]);

  inter2 = gridfit(X(:), Y(:), signals(:), p, t, 'smooth', [25 500]);
  flow_full = interp2(X(:,:,1), Y(:,:,1), inter2, [1:66], t(goods));
  flow_full = padarray(flow_full, [10 1], 0);

  flow_small = interp2(X(:,:,1), Y(:,:,1), inter2, [1:66], [range_t(1):15:range_t(end)].');
  flow_small = padarray(flow_small, [1 1], 0);

  imagesc(flow_full);
  set(gca, 'XTickLabel', pos(get(gca, 'XTick')));
  colorbar;
  print('-dpng', ['PNG/flow_left_smooth.png']);

  signals = max(signals_full(:,(end/2)+1:end,:), - signals_full(:,end/2:-1:1,:));

  valids = (sum(~isnan(signals), 3) <= 2);
  signals(repmat(valids, [1 1 nfiles])) = 0;

  signals = padarray(signals, [10 3 0], 0);
  p = pos((end/2)+1:end);
  t = [0:size(signals, 1)-1].';
  p = [[-2:0].'; p; [67:69].'];

  [X,Y] = meshgrid(p, t);
  X = repmat(X, [1 1 nfiles]);
  Y = repmat(Y, [1 1 nfiles]);

  goods = any(~valids, 2);
  range_t = t(goods);
  range_t = range_t([1 end]);

  inter2 = gridfit(X(:), Y(:), signals(:), p, t, 'smooth', [25 500]);
  flow_full = interp2(X(:,:,1), Y(:,:,1), inter2, [1:66], t(goods));
  flow_full = padarray(flow_full, [10 1], 0);

  flow_small = interp2(X(:,:,1), Y(:,:,1), inter2, [1:66], [range_t(1):15:range_t(end)].');
  flow_small = padarray(flow_small, [1 1], 0);

  imagesc(flow_full);
  set(gca, 'XTickLabel', pos(get(gca, 'XTick')));
  colorbar;
  print('-dpng', ['PNG/flow_max_smoothed.png']);


  keyboard

  flow = fliplr(flow_full);
  flow2 = gaussian_mex(flow, 0.67);
  flow3 = gaussian_mex(flow2, 0.67);

  flow_step = 1;
  save('cyto_flow', 'flow', 'flow2', 'flow3', 'flow_step');

  flow = fliplr(flow_small);
  flow2 = gaussian_mex(flow, 0.67);
  flow3 = gaussian_mex(flow2, 0.67);
  flow_step = 15;
  save('cyto_flow_small', 'flow', 'flow2', 'flow3', 'flow_step');

  %set(gca, 'XTickLabel', res{end,3}(get(gca, 'XTick')));
  %  colorbar;
  %print('-dpng', '-r450', ['PNG/flow_avg_2.png']);
  
  %[m,s] = (mymean(signals(end/2+1:end,:,:)-signals(end/2:-1:1,:,:), 3));
  %m(s==0) = NaN;
  %imagesc(m.');
  %set(gca, 'XTickLabel', res{end,3}(get(gca, 'XTick')+end/2));
  %  colorbar;
  %print('-dpng', '-r450', ['PNG/flow_avg_half_2.png']);
  %subplot(2,2,4);

  return;

  keyboard

  flow_files = dir('NMY2Flow_NMY2_*.mat');

  x_bin = 4;
  y_bin = 20;

  XI = [];
  YI = [];
  ZI = [];
  for i=1:length(flow_files)
    f = load(flow_files(i).name);

    t_XI = f.XI;
    t_YI = f.YI;
    t_ZI = f.ZI;

    t_ZI(t_XI == 0) = NaN;
    t_ZI(t_XI == min(t_XI(:))) = NaN;

    x = unique(t_XI);
    y = unique(t_YI);
    figure;imagesc(-t_ZI);
    set(gca, 'XTickLabel', round(x(get(gca, 'XTick'))), 'YTickLabel', round(y(get(gca, 'YTick'))));
    print('-dpng', '-r450', ['PNG/flow_nate_' num2str(i) '.png']);

    goods = isfinite(t_ZI);
    t_XI = t_XI(goods);
    t_YI = t_YI(goods);
    t_ZI = t_ZI(goods);

    XI = [XI; t_XI];
    YI = [YI; t_YI];
    ZI = [ZI; t_ZI];
  end

  x = [-67:x_bin:0 0];
  y = [min(YI):y_bin:max(YI) max(YI)];

  x(1) = x(1) - 1e-5;
  x(end) = x(end) + 1e-5;

  y(1) = y(1) - 10;
  y(end) = y(end) + 10;

  [map, edges, centers, pos] = histcn([XI YI], x, y);
  [coords, junk, index] = unique(pos, 'rows');

  vals = NaN(size(coords, 1), 1);
  full_vals = NaN(size(map));
  for i=1:size(coords, 1)
    currents = (index==i);
    full_vals(coords(i,1), coords(i,2)) = mymean(ZI(currents));
  end

  [X, Y] = meshgrid([centers{1}(1)-x_bin centers{1} centers{1}(end)+x_bin], ...
                    [centers{2}(1)-y_bin centers{2} centers{2}(end)+y_bin]);

  x = unique(X);
  y = unique(Y);

  full_vals = padarray(full_vals, [1 1], 0).';
  goods = isfinite(full_vals);
  Z = griddata(X(goods), Y(goods), full_vals(goods), X, Y);
  figure;imagesc(-Z);
  set(gca, 'XTickLabel', round(x(get(gca, 'XTick'))), 'YTickLabel', round(y(get(gca, 'YTick'))));
  print('-dpng', '-r450', ['PNG/flow_nate_avg-' num2str(x_bin) '.png']);

  full_vals(Y >= 0) = NaN;
  full_vals(Y == max(Y(:))) = 0;
  goods = isfinite(full_vals);
  Z_trunk = griddata(X(goods), Y(goods), full_vals(goods), X, Y);
  figure;imagesc(-Z_trunk);
  set(gca, 'XTickLabel', round(x(get(gca, 'XTick'))), 'YTickLabel', round(y(get(gca, 'YTick'))));
  print('-dpng', '-r450', ['PNG/flow_nate_cut-' num2str(x_bin) '.png']);

  return;
end
