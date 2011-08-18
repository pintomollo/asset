function [bests, indxs] =  dp_score_2d(datas, prev_datas, prev_dist, prev_dir, initiation, params)

  %keyboard

  half = params.nhood;
  alpha = params.alpha;
  beta = params.beta;
  gamma = params.gamma;
  delta = params.delta;

  if (numel(half) == 1)
    half = [half half];
  end
  if (numel(alpha) == 1)
    alpha = [alpha alpha];
  end
  if (numel(beta) == 1)
    beta = [beta beta];
  end
  if (numel(gamma) == 1)
    gamma = [gamma gamma];
  end
  if (numel(delta) == 1)
    delta = [delta delta];
  end

  if (sum(delta) == 0)
    delta = delta + 1;
  end
  delta = delta / sum(delta);

  %keyboard

  prohibit = params.prohibit;
  new_prctile = params.spawn_percentile * 100;
  max_nhood = max(half);

  %half = floor(nhood/2);
  %check = [-half:half];
  %nhood = (check)^2;

  [tmpx, tmpy] = meshgrid([-half(1):half(1)], [-half(2):half(2)]);

  %dists = sqrt(bsxfun(@plus, ([-half(1):half(1)].').^2, [-half(2):half(2)].^2));
  window_size = size(tmpx.');
  nhood = numel(tmpx);
  window_indx = [1:nhood].';

  dists = abs([tmpx(:) tmpy(:)]);

  %if (half~=0)
  max_dist = max(dists, [], 1);
  max_dist(max_dist == 0) = 1;
  dists = bsxfun(@rdivide, dists, max_dist);
  weighted_dists = sum(bsxfun(@times, dists, alpha .* beta .* gamma .* delta), 2);

  %  dists = dists/max(dists);
  %end

  weighted_datas = datas * sum((1-alpha) .* delta);
  weighted_initiation = initiation * sum(alpha .* delta);

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

  real_indx = curr_indx(1:window_size(1),1:window_size(2));
  real_indx = real_indx(:);
  back_indx = spalloc(npts, 1, nhood);
  back_indx(real_indx) = window_indx;

  %curr_dir = curr_indx - 

  res = zeros([nhood,size_prob]);
  if (isempty(prev_dist))
    prev_dist = zeros(size(initiation));
  end
      prev_dist = mirror_matrix(prev_dist, half);

  %if (~isempty(prev_dir))
  %  new_segments = (prev_dir == 0);
  %  prev_dir(new_segments) = curr_indx(new_segments);

  %  tmp = bsxfun(@minus, prev_dir(:) - curr_indx(:), [0 npts -npts]);
  %  [~, indxs] = min(abs(tmp), [], 2);
  %  prev_dir = tmp(sub2ind([npts, 3], curr_indx(:), indxs)) + real_indx(ceil(nhood/2));
  %  prev_dir = reshape(prev_dir, size_prob);
  %  prev_dir = back_indx(prev_dir);
    
    
    %[dir_i, dir_j] = ind2sub(npts, prev_dir);

    %[~, indx] = min(cat(3, abs(bsxfun(@minus, dir_i, [1:npts(1)].')), abs(bsxfun(@minus, dir_j, [1:npts(2)]))),[],3);

    %tmp = [prev_dir - [1:npts]; prev_dir - [1:npts] - npts];
    %[junk, indxs] = min(abs(tmp));

    %prev_dir = tmp(sub2ind([2,npts],indxs,[1:npts]));
  %end

  %uppers = mirror_matrix(uppers, half);
  %prev_val = [Inf*ones(1,half) values Inf*ones(1,half)];
  %prev_dir = [Inf*ones(1,half) prev_dir Inf*ones(1,half)];
  if (~isempty(prev_datas))
    %prev_val = mirror_matrix(values, nhood);
    %prev_val = values(1, [end-max_nhood+1:end 1:end 1:max_nhood]);
    prev_datas = mirror_matrix(prev_datas, half);

    new_segments = (prev_dir == 0);
    prev_dir(new_segments) = curr_indx(new_segments);

    tmp = bsxfun(@minus, prev_dir(:) - curr_indx(:), [0 npts -npts]);
    [~, indxs] = min(abs(tmp), [], 2);
    prev_dir = tmp(sub2ind([npts, 3], curr_indx(:), indxs)) + real_indx(ceil(nhood/2));
    prev_dir = reshape(prev_dir, size_prob);
    prev_dir = back_indx(prev_dir);
    prev_dir = mirror_matrix(prev_dir, half);

  %if (nargout > 2)
    %  do_probs = true;
  %
  %    emission = zeros(nhood,npts);
  %  else
  %    do_probs = false;
  %  end
    %size_prob = size_prob(1);
    tmp_inf = Inf(1, max(size_prob));
    max_dir = 2*max(dists, [], 1);
    weighted_max = (alpha .* beta .* (1-gamma) .* delta) ./ max_dir;
    [tmp_i, tmp_j] = ind2sub(window_size, window_indx);
    curr_dir = [tmp_i tmp_j];
    crossing = zeros(size_prob);

    for i=1:nhood


      switch prohibit
        case 'diag'

          diag_indx = (tmp_j(i) - tmp_i(i));
          diag_switch = sign(diag_indx);

          if (diag_switch > 0)
            diag_indx = [-diag_indx+1 : 0];
          elseif (diag_switch < 0)
            diag_indx = [0:-diag_indx-1];
          else
            diag_indx = [];
          end

          crossing = zeros(size_prob);
          for j=diag_indx
            crossing = crossing + diag(tmp_inf(1:end-abs(j)), j);
          end
      end
      %disp([num2str(prev_val(i))  ' '  num2str(new_values(1))  ' '  num2str(datas(1))  ' '  num2str(prev_dist(i)) ' ' num2str(prev_dir(i))])

      if(isempty(prev_datas))
        weighted_smooth = 0;
      else
        weighted_dir = dir_dist(curr_dir(i, :),prev_dir(tmp_i(i):tmp_i(i)+size_prob(1)-1, tmp_j(i):tmp_j(i)+size_prob(2)-1), window_size, weighted_max);

        weighted_val = (abs(datas - prev_datas(tmp_i(i):tmp_i(i)+size_prob(1)-1, tmp_j(i):tmp_j(i)+size_prob(2)-1)) / 2) * sum(alpha .* (1-beta) .* delta);


        weighted_smooth = weighted_dists(i) + weighted_dir + weighted_val;
      end

      res(i,:,:) = weighted_smooth + weighted_datas + ...
            prev_dist(tmp_i(i):tmp_i(i)+size_prob(1)-1,tmp_j(i):tmp_j(i)+size_prob(2)-1) + ...
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
      %keyboard

      %weighted_initiation = sum(bsxfun(@times, initiation, reshape(alpha .* delta, [1 1 2])), 3);
      starts = weighted_initiation + weighted_datas + prctile(prev_dist(~isinf(prev_dist)), new_prctile) + crossing;
      [bests, tmp] = min(cat(3,bests,starts), [], 3);
    end
  else
    bests = weighted_initiation + weighted_datas;
    indxs = curr_indx;
    tmp = 2*ones(size(indxs));
  end

  if (all(isinf(bests)))
    %disp('Warning : No transition is valid !!');
    indxs = curr_indx;
  end

  indxs(indxs < 1) = indxs(indxs < 1) + npts;
  indxs(indxs > npts) = indxs(indxs > npts) - npts;
  if (new_prctile >= 0)
    indxs(tmp == 2) = 0;
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

  npts = size(matrix);
  new_matrix = NaN(npts + 2*nhood);

  new_matrix([nhood(1)+1:end-nhood(1)],[nhood(2)+1:end-nhood(2)]) = matrix;
  new_matrix(:,[1:nhood(2)]) = matrix([end-nhood(1)+1:end 1:end 1:nhood(1)], [end-nhood(2)+1:end]);
  new_matrix(:,[end-nhood(2)+1:end]) = matrix([end-nhood(1)+1:end 1:end 1:nhood(1)], [1:nhood(2)]);
  new_matrix([1:nhood(1)],[nhood(2)+1:end-nhood(2)]) = matrix([end-nhood(1)+1:end], :);
  new_matrix([end-nhood(1)+1:end],[nhood(2)+1:end-nhood(2)]) = matrix([1:nhood(1)], :);

  return;
end

function vals = dir_dist(curr_dir, dirs, window_size, weight)

  center_pix = ceil(prod(window_size)/2);

  [dir_i, dir_j] = ind2sub(window_size, dirs);
  %vals = sqrt((curr_dir(1) - dir_i).^2 + (curr_dir(2) - dir_j).^2);
  vals = abs(curr_dir(1) - dir_i)*weight(1) + abs(curr_dir(2) - dir_j)*weight(2);

  % New paths have no direction
  vals(dirs == center_pix) = 0;

  return;
end
