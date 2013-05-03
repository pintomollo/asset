function [flows, flows_small, pnm] = average_flow(signals, alignment, proj_bin_size)

  if (ischar(signals))
    files = dir(signals);
    nfiles = length(files);
    signals = cell(nfiles, 1);
    pnm = NaN(nfiles, 1);

    for i=1:nfiles
      load(files(i).name);
      time = get_manual_timing(mymovie, opts);
      pnm(i) = time(2);
      tmp = display_flow(mymovie, opts);
      signals{i} = tmp.';
    end
  elseif (~iscell(signals))
    imgs = signals;
    nimgs = size(signals, 3);
    signals = cell(nimgs, 1);

    for i=1:nimgs
      img = imgs(:,:,i);
      
      bads = isnan(img);
      cols = all(bads, 1);
      rows = all(bads, 2);

      img = img(~rows, :);
      img = img(:, ~cols);

      signals{i} = img;
    end
    pnm = zeros(nimgs, 1);
  else
    pnm = zeros(length(signals), 1);
  end

  [signals, shift] = stack_images(signals, alignment);
  half = size(signals, 2) / 2;

  avg = nanmean(signals, 3);
  clim = [min(avg(:)) max(avg(:))];
  for i=1:nfiles
    figure;imagesc(signals(:,:,i), clim);
  end

  figure;imagesc(avg);

  pos = [-half:half]*proj_bin_size;
  pos = (pos(1:end-1) + pos(2:end)) / 2;

  %signals = (signals(:,half+1:end,:) - signals(:,half:-1:1,:)) / 2;
  signals = cat(3,signals(:,half+1:end,:), -signals(:,half:-1:1,:));
  all_signals = signals;

  nimgs = size(signals, 3);

  counts = nanmean(isfinite(signals), 3);
  bads = (counts < 0.5);
  pnm = pnm + shift;

  %signals(repmat(bads, [1 1 nimgs])) = 0;
  signals(~isfinite(signals)) = 0;

  signals = padarray(signals, [100 0 0], 0);
  counts = padarray(counts, [100 0 0], 0);
  t = [0:size(signals, 1)-1].';
  p = pos(half+1:end).';
  %left = [-2:p(1)-proj_bin_size/2];
  %right = [p(end)+proj_bin_size/2 : 69];
  left = [-10:p(1)-proj_bin_size/2];
  right = [p(end)+proj_bin_size/2 : 85];

  p = [left(:); p; right(:)];
  signals = padarray(signals, [0 length(left) 0], 0, 'pre');
  signals = padarray(signals, [0 length(right) 0], 0, 'post');
  counts = padarray(counts, [0 length(left) 0], 0, 'pre');
  counts = padarray(counts, [0 length(right) 0], 0, 'post');

  [X,Y] = meshgrid(p, t);
  X = repmat(X, [1 1 nimgs]);
  Y = repmat(Y, [1 1 nimgs]);

  inter2 = gridfit(X(:), Y(:), signals(:), p, t, 'smooth', [25 500]);

  %% Use counts to make a mask with (>0.5) and smooth it linearly or exponentially,
  %% weight the speed to put them to 0, keep the moving inner part

  counts = (counts > 0.5);
  counts = double(bwdist(counts));
  weights = exp(-counts/(2*5^2));

  inter2 = inter2 .* weights;

  %goods = any(inter2 ~= 0, 2);
  goods = (abs(mean(inter2, 2)) > max(inter2(:))/1000);
  range_t = t(goods);
  range_t = range_t([1 end]);

  flow_full = interp2(X(:,:,1), Y(:,:,1), inter2, [1:66], t(goods));
  flow_full = padarray(flow_full, [10 1], 0);

  flow_small = interp2(X(:,:,1), Y(:,:,1), inter2, [1:66], [range_t(1):15:range_t(end)].');
  flow_small = padarray(flow_small, [1 1], 0);

  flow = fliplr(flow_full);
  flow2 = gaussian_mex(flow, 0.67);
  flow3 = gaussian_mex(flow2, 0.67);

  flows = cat(3, flow, flow2, flow3);

  flow_step = 1;

  flow = fliplr(flow_small);
  flow2 = gaussian_mex(flow, 0.67);
  flow3 = gaussian_mex(flow2, 0.67);
  flow_step = 15;

  flows_small = cat(3, flow, flow2, flow3);

  return;
end
