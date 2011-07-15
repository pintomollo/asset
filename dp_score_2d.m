function [bests, indxs] =  dp_score_2d(values, new_values, datas, prev_dist, prev_dir, initiation, params)

  %keyboard

  half = params.nhood;
  alpha = params.alpha;
  beta = params.beta;
  gamma = params.gamma;

  prohibit = params.prohibit;
  new_prctile = params.spawn_percentile;

  %half = floor(nhood/2);
  check = [-half:half];
  %nhood = (check)^2;

  dists = sqrt(bsxfun(@plus, (check.').^2, check.^2));
  window_size = size(dists);
  dists = dists(:);
  nhood = numel(dists);
  window_indx = [1:nhood].';

  if (half~=0)
    dists = dists/max(dists);
  end

  %%%% TO DO
  % Identify the corner of useful data (find(~isinf(datas),1))
  % Take out central part
  % Run normally
  % re-insert results

  size_prob = size(datas);
  npts = numel(datas);
  curr_indx = reshape([1:npts], size_prob);

  %[indxi, indxj] = ind2sub(size_prob,curr_indx);
  %uppers = (indxi <= indxj);

  real_indx = curr_indx(1:window_size(1),1:window_size(1));
  real_indx = real_indx(:);
  back_indx = spalloc(npts, 1, nhood);
  back_indx(real_indx) = window_indx;

  %curr_dir = curr_indx - 

  res = zeros([nhood,size_prob]);
  new_segments = (prev_dir == 0);
  prev_dir(new_segments) = curr_indx(new_segments);

  if (~isempty(prev_dir))
    tmp = bsxfun(@minus, prev_dir(:) - curr_indx(:), [0 npts -npts]);
    [~, indxs] = min(abs(tmp), [], 2);
    prev_dir = tmp(sub2ind([npts, 3], curr_indx(:), indxs)) + real_indx(ceil(nhood/2));
    prev_dir = reshape(prev_dir, size_prob);
    prev_dir = back_indx(prev_dir);
    
    
    %[dir_i, dir_j] = ind2sub(npts, prev_dir);

    %[~, indx] = min(cat(3, abs(bsxfun(@minus, dir_i, [1:npts(1)].')), abs(bsxfun(@minus, dir_j, [1:npts(2)]))),[],3);

    %tmp = [prev_dir - [1:npts]; prev_dir - [1:npts] - npts];
    %[junk, indxs] = min(abs(tmp));

    %prev_dir = tmp(sub2ind([2,npts],indxs,[1:npts]));
  end

  prev_dist = mirror_matrix(prev_dist, half);
  %uppers = mirror_matrix(uppers, half);
  %prev_val = [Inf*ones(1,half) values Inf*ones(1,half)];
  %prev_dir = [Inf*ones(1,half) prev_dir Inf*ones(1,half)];
  if (~isempty(values))
    %prev_val = mirror_matrix(values, nhood);
    prev_val = values(1, [end-half+1:end 1:end 1:half]);
    prev_dir = mirror_matrix(prev_dir, half);
  end

  %if (nargout > 2)
  %  do_probs = true;
%
%    emission = zeros(nhood,npts);
%  else
%    do_probs = false;
%  end
  size_prob = size_prob(1);
  tmp_inf = Inf(1, size_prob);
  max_dir = 2*max(dists);
  [tmp_i, tmp_j] = ind2sub(window_size, window_indx);
  curr_dir = [tmp_i tmp_j];
  crossing = zeros(size_prob);

  for i=1:nhood
    diag_indx = (tmp_j(i) - tmp_i(i));
    diag_switch = sign(diag_indx);

    if (diag_switch > 0)
      diag_indx = [-diag_indx+1 : 0];
    elseif (diag_switch < 0)
      diag_indx = [0:-diag_indx-1];
    else
      diag_indx = [];
    end

    switch prohibit
      case 'diag'
        crossing = zeros(size_prob);
        for j=diag_indx
          crossing = crossing + diag(tmp_inf(1:end-abs(j)), j);
        end
    end
    %disp([num2str(prev_val(i))  ' '  num2str(new_values(1))  ' '  num2str(datas(1))  ' '  num2str(prev_dist(i)) ' ' num2str(prev_dir(i))])

    if(isempty(values))
      smooth = 0;
    else
      smooth = (dists(i)*gamma + ...
        (dir_dist(curr_dir(i, :),prev_dir(tmp_i(i):tmp_i(i)+size_prob-1, tmp_j(i):tmp_j(i)+size_prob-1), window_size)/max_dir)*(1-gamma))*beta + ...
        (bsxfun(@plus, abs(new_values - prev_val(tmp_i(i):tmp_i(i)+size_prob-1)).', ...
        abs(new_values - prev_val(tmp_j(i):tmp_j(i)+size_prob-1))) / 2)*(1-beta);
    end

    res(i,:,:) = alpha*smooth + (1-alpha)*datas + ...
          prev_dist(tmp_i(i):tmp_i(i)+size_prob-1,tmp_j(i):tmp_j(i)+size_prob-1) + ...
          crossing;
    %      Inf*xor(uppers(half+1:half+size_prob, half+1:half+size_prob), uppers(tmp_i(i):tmp_i(i)+size_prob-1,tmp_j(i):tmp_j(i)+size_prob-1));


%    if (do_probs)
%      if(isempty(values))
%        emission(i,:) = (1-alpha)*datas;
%      else
%        emission(i,:) = alpha * ((abs(check(i)-prev_dir(i:i+npts-1))/(2*half))*(1-gamma)*beta + abs(new_values-prev_val(i:i+npts-1))*(1-beta)) + (1-alpha)*datas;
%      end
%    end
  end
  res(isnan(res)) = Inf;

  [bests, indxs] = min(res,[],1);
  bests = squeeze(bests);
  indxs = squeeze(indxs);
  direction = real_indx(indxs) - real_indx(ceil(nhood/2));
  %direction = indxs - ceil(nhood/2);
  indxs = curr_indx + direction;

  if (new_prctile >= 0)
    starts = alpha*initiation + (1-alpha)*datas + prctile(prev_dist(:), new_prctile);
    [bests, tmp] = min(cat(3,starts,bests), [], 3);
  end

  if (all(isinf(bests)))
    %disp('Warning : No transition is valid !!');
    indxs = curr_indx;
  end

  indxs(indxs < 1) = indxs(indxs < 1) + npts;
  indxs(indxs > npts) = indxs(indxs > npts) - npts;
  if (new_prctile >= 0)
    indxs(tmp == 1) = 0;
  end

%  if (do_probs)
%
%    emission(isnan(emission)) = Inf;
%    emission = exp(-emission);
%
%    if (nargout > 3)
%      if(isempty(values))
%        trans = 1;
%      else
%        trans = exp(-dists*alpha*beta*gamma);
%      end
%    end
%  end

  return;
end

function new_matrix = mirror_matrix(matrix, nhood)

  npts = size(matrix, 1);
  new_matrix = NaN(npts + 2*nhood);

  new_matrix([nhood+1:end-nhood],[nhood+1:end-nhood]) = matrix;
  new_matrix(:,[1:nhood]) = matrix([end-nhood+1:end 1:end 1:nhood], [end-nhood+1:end]);
  new_matrix(:,[end-nhood+1:end]) = matrix([end-nhood+1:end 1:end 1:nhood], [1:nhood]);
  new_matrix([1:nhood],[nhood+1:end-nhood]) = matrix([end-nhood+1:end], :);
  new_matrix([end-nhood+1:end],[nhood+1:end-nhood]) = matrix([1:nhood], :);

  return;
end

function vals = dir_dist(curr_dir, dirs, window_size)

  [dir_i, dir_j] = ind2sub(window_size, dirs);
  vals = sqrt((curr_dir(1) - dir_i).^2 + (curr_dir(2) - dir_j).^2);

  return;
end
