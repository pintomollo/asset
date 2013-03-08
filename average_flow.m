function [flows, flows_small] = average_flow(signals, alignment, proj_bin_size)

  if (~iscell(signals))
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
  end

  signals = stack_images(signals, alignment);
  half = size(signals, 2) / 2;

  pos = [-half:half]*proj_bin_size;
  pos = (pos(1:end-1) + pos(2:end)) / 2;

  %signals = (signals(:,half+1:end,:) - signals(:,half:-1:1,:)) / 2;
  signals = cat(3,signals(:,half+1:end,:), -signals(:,half:-1:1,:));

  nimgs = size(signals, 3);

  counts = nanmean(isfinite(signals), 3);
  bads = (counts < 2/3);
  %signals(repmat(bads, [1 1 nimgs])) = 0;
  signals(~isfinite(signals)) = 0;

  signals = padarray(signals, [10 0 0], 0);
  t = [0:size(signals, 1)-1].';
  p = pos(half+1:end).';
  left = [-2:p(1)-proj_bin_size/2];
  right = [p(end)+proj_bin_size/2 : 69];
  p = [left(:); p; right(:)];
  signals = padarray(signals, [0 length(left) 0], 0, 'pre');
  signals = padarray(signals, [0 length(right) 0], 0, 'post');

  [X,Y] = meshgrid(p, t);
  X = repmat(X, [1 1 nimgs]);
  Y = repmat(Y, [1 1 nimgs]);

  goods = any(~bads, 2);
  range_t = t(goods);
  range_t = range_t([1 end]);

  inter2 = gridfit(X(:), Y(:), signals(:), p, t, 'smooth', [25 500]);
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
