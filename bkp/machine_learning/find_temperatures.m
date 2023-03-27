function [beta, gamma, map] = find_temperatures(trans, emission, opts)

  sparse_thresh = 1e-3;

  [nindx, nhood, npts] = size(emission);

  aim_dist = opts.aim_transitions * npts / 2;
  aim_stds = opts.aim_emissions * npts / 2;
  thresh = opts.thresh;
  step_thresh = opts.step_thresh;

  beta = 1;
  beta_range = [0 Inf];
  gamma = 1;
  gamma_range = [0 Inf];

  inner_trans = trans(2:end-1);
  half_dist = floor((npts + 2) / 2);
  dist = [([1:half_dist] - 1) (npts + 1 - [half_dist+1:npts])];

  prev_dist = Inf;
  exp_dist = expected_deviation(inner_trans, gamma, npts, nindx, nhood, dist);
  while (abs(exp_dist - aim_dist) > thresh && abs(prev_dist - exp_dist) > step_thresh)
    if (exp_dist < aim_dist)
      gamma_range(2) = gamma;
      gamma = mean(gamma_range);
    else
      gamma_range(1) = gamma;
      if (isinf(gamma_range(2)))
        gamma = 10 * gamma;
      else
        gamma = mean(gamma_range);
      end
    end
    prev_dist = exp_dist;
    exp_dist = expected_deviation(inner_trans, gamma, npts, nindx, nhood, dist);
  end

  prev_stds = Inf;
  [map, mean_stds] = posterior_decoding(emission, trans, beta, gamma);
  while (abs(mean_stds - aim_stds) > thresh && abs(prev_stds - mean_stds) > step_thresh)
    if (mean_stds < aim_stds)
      beta_range(2) = beta;
      beta = mean(beta_range);
    else
      beta_range(1) = beta;
      if (isinf(beta_range(2)))
        beta = 10 * beta;
      else
        beta = mean(beta_range);
      end
    end
    prev_stds = mean_stds;
    [map, mean_stds] = posterior_decoding(emission, trans, beta, gamma);
  end

  map(map < sparse_thresh) = 0;
  map = map ./ repmat(sum(map,2),1,npts);
  %map = sparse(map);
  %values = [[0:0.1:1] [1:1:10] [10:10:100] [100:100:1000]]
  %values = [1:10:1000];
  %entr = [];
  %stds = [];

  %for i = values
  %  i
  %  [entr(end+1), stds(end+1)] = posterior_decoding([], emission, trans, i, gamma);
  %end

  %keyboard
  %corrs
  %msd = sum(corrs)

  %implot(map)

  return;
end

function exp_dist = expected_deviation(trans, beta, npts, nindx, nhood, dist)

  even = (mod(nhood,2) == 0);
  half = floor(nhood / 2);

  if (beta == 0)
    trans = ones(size(trans));
  else
    trans = log(trans) * beta;
    trans = exp(trans - max(trans));
    trans = trans / sum(trans);
  end

  %trans = [0 0 0 1 0 0 0]

  map = zeros(nindx,npts);
  distrib = zeros(1,npts);
  distrib(1,1) = 1;
  %corrs = zeros(nindx,1);

  for i=1:nindx
    tmp = cconv(distrib, trans, npts);
    distrib = tmp([half+1:end 1:half]);
    if (even && (mod(i,2)==0))
      distrib = distrib([end 1:end-1]);
    end
    map(i,:) = distrib;
    %corrs(i) = corr(map(1,:).',map(i+1,:).');
  %end
  %for i=1:nindx
  %  corrs(i) = abs(corr(map(2,:).',map(i+1,:).')) * (nindx - i - 1);
  end
  exp_dist = sum(dist .* distrib);
  %figure;implot(map)

  return;
end
