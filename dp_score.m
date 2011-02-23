function [bests, indxs, emission, trans] =  dp_score(values, candidates, datas, prev_dist, prev_dir, params)

  nhood = params.nhood;
  alpha = params.alpha;
  beta = params.beta;
  gamma = params.gamma;

  half = floor(nhood/2);
  check = [-half:half];

  if (half==0)
    dists = abs(check);
  else
    dists = abs(check)/half;
  end

  npts = length(candidates);
  res = zeros(nhood,npts);

  if (~isempty(prev_dir))
    tmp = [prev_dir - [1:npts]; prev_dir - [1:npts] - npts];
    [junk, indxs] = min(abs(tmp));
    prev_dir = tmp(sub2ind([2,npts],indxs,[1:npts]));
  end

  prev_dist = prev_dist(1, [end-half+1:end 1:end 1:half]);
  %prev_dist = [Inf*ones(1,half) prev_dist Inf*ones(1,half)];
  %prev_val = [Inf*ones(1,half) values Inf*ones(1,half)];
  %prev_dir = [Inf*ones(1,half) prev_dir Inf*ones(1,half)];
  if (~isempty(values))
    prev_val = values(1, [end-half+1:end 1:end 1:half]);
    prev_dir = prev_dir(1, [end-half+1:end 1:end 1:half]);
  end

  if (nargout > 2)
    do_probs = true;

    emission = zeros(nhood,npts);
  else
    do_probs = false;
  end

  for i=1:nhood
    %disp([num2str(prev_val(i))  ' '  num2str(candidates(1))  ' '  num2str(datas(1))  ' '  num2str(prev_dist(i)) ' ' num2str(prev_dir(i))])

    if(isempty(values))
      smooth = 0;
    else
      smooth = (dists(i)*gamma + (abs(check(i)-prev_dir(i:i+npts-1))/(2*half))*(1-gamma))*beta + abs(candidates-prev_val(i:i+npts-1))*(1-beta);
    end

    res(i,:) = alpha*smooth + (1-alpha)*datas + prev_dist(i:i+npts-1);

    res(i,isnan(res(i,:))) = Inf;

    if (do_probs)
      if(isempty(values))
        emission(i,:) = (1-alpha)*datas;
      else
        emission(i,:) = alpha * ((abs(check(i)-prev_dir(i:i+npts-1))/(2*half))*(1-gamma)*beta + abs(candidates-prev_val(i:i+npts-1))*(1-beta)) + (1-alpha)*datas;
      end
    end
  end

  [bests, indxs] = min(res,[],1);
  direction = indxs - half - 1;
  indxs = [1:npts] + direction;

  if (all(isinf(bests)))
    %disp('Warning : No transition is valid !!');
    indxs = [1:npts];
  end

  indxs(indxs < 1) = indxs(indxs < 1) + npts;
  indxs(indxs > npts) = indxs(indxs > npts) - npts;

  if (do_probs)

    emission(isnan(emission)) = Inf;
    emission = exp(-emission);

    if (nargout > 3)
      if(isempty(values))
        trans = 1;
      else
        trans = exp(-dists*alpha*beta*gamma);
      end
    end
  end

  return;
end
