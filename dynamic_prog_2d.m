function [path] = dynamic_prog_2d(img, params, weight, weight_params, opts)

  [h,w] = size(img);

  nhood = params.nhood;
  half = floor(nhood/2);

  indxs = 1:h;

  nsteps = length(indxs);

  init = params.init; 
  
  path = zeros(h,1);
  if (isempty(weight_params))
    wimg = weight(img);
  else
    wimg = weight(img, weight_params);
  end

  size_weight = size(wimg);
  size_problem = size_weight(2:end);

  dist = zeros([nsteps size_problem]);
  map = zeros([nsteps size_problem]);
  %if (do_probs)
  %  emission = zeros(nsteps, nhood, w);
  %  transitions = zeros(1,nhood+2);
  %end

  if (isempty(init))
    init = squeeze(dist(1,:,:));
  elseif (numel(init) < numel(dist(1,:,:)))
    tmp = init;
    init = Inf(size(squeeze(dist(1,:,:))));
    init(tmp) = 0;
  end
  %nstart = sum(~isinf(init(:)));

  prev_step = img(1,:);
  %if (do_probs)
  %  [dist(1,:), map(1,:), emission(1,:,:), transitions(1,1)] = dp_score_mex([], prev_line, wimg(1,:), init, [], params);
  %else
    [dist(1,:,:), map(1,:,:)] = dp_score_2d([], prev_step, squeeze(wimg(1,:,:)), init, [], params);
  %end

  for j=2:nsteps
    indx = indxs(j)

    step = img(indx, :);

    %if (do_probs)
    %  if (j==2)
    %    [dist(j,:), map(j,:), emission(j,:,:), transitions(1,2:end-1)] = dp_score_mex(prev_line, line, wimg(indx,:), dist(j-1,:), map(j-1,:), params);
    %  else
    %    [dist(j,:), map(j,:), emission(j,:,:)] = dp_score_mex(prev_line, line, wimg(indx,:), dist(j-1,:), map(j-1,:), params);
    %  end
    %else
      [dist(j,:,:), map(j,:,:)] = dp_score_2d(prev_step, step, squeeze(wimg(indx,:,:)), squeeze(dist(j-1,:,:)), squeeze(map(j-1,:,:)), params);
    %end

    if (all(all(isinf(dist(j,:,:)))))
      dist(j,:,:) = dist(j-1,:,:);
      map(j,:,:) = reshape([1:numel(step)], size_problem);
    end

    prev_step = step;
  end

  %keyboard

  %if (do_probs)
  %  [junk, junk, junk, transitions(1,end)] = dp_score_mex([], prev_line, wimg(indxs(end),:), init, [], params);
  %end

  tmp = zeros(nsteps, 1);

  %if (isfield(params, 'final') & ~isempty(params.final))
  %  subdist = dist(end, params.final);
  %  tmp(end,1) = params.final(find(subdist == min(subdist), 1));
  %else
    tmp(end,1) = find(dist(end,:,:) == min(min(dist(end,:,:))),1);
  %end

  for j=1:nsteps-1
    [tmpi, tmpj] = ind2sub(size_problem,tmp(end-j+1,1));
    tmp(end-j,1) = map(end-j+1,tmpi, tmpj);
  end

%  switch opts.dp_method
%    case 'backward'
%      path = tmp(end:-1:h,1);
%      if (do_probs)
%        emission = emission(end:-1:h, : , :);
%      end
%    case 'double'
%      middle = floor(h/2);
%      path = tmp([h+1:h+middle middle+1:h],1);
%
%      if (do_probs)
%        emission = emission([h+1:middle+h middle+1:h], :, :);
%      end
%    otherwise
    [tmpi, tmpj] = ind2sub(size_problem,tmp);
      path = [tmpi(:) tmpj(:)];
%  end

  %if (opts.force_circularity && abs(path(1,1)-path(end,1))>(nhood/2))
  %  
  %  check = [-half:half];
%
%    if (nstart ~= 1)

%      if (nstart < w & nstart > 0)
%        dcheck = floor(nstart / 2) - 1;
%        check = [-dcheck:dcheck];
%
%        if (isempty(check))
%          check = 0;
%        end
%      end
%      check = check + path(end,1);
%      check = check(check > 0 & check <= w);

%      init = Inf(1,w);
%      init(1, check) = 0;
%
%      params.init = init;

%      if (do_probs)
%        [path, emission, transitions] = dynamic_programming(img, params, weight, weight_params, opts);
%      else
%        [path] = dynamic_programming(img, params, weight, weight_params, opts);
%      end
%    else
%      check = check + path(1,1);
%      check = check(check > 0 & check <= w);
%
%      [junk, mins] = min(dist(1, check));
%
%      path(end,1) = check(mins);
%
%      for j=1:h-1
%        new_path = map(end-j+1,path(end-j+1,1));
%        if (path(end-j,1)==new_path)
%          break;
%        end
%        path(end-j,1) = new_path;
%      end
%    end
%  end

%  if (retrieve_maps)
%    switch opts.dp_method
%      case 'backward'
%        emission = map(end:-1:h,:);
%        transitions = dist(end:-1:h,:);
%      case 'double'
%        middle = floor(h/2);
%        emission = map([h+1:h+middle middle+1:h],:);
%        transitions = dist([h+1:h+middle middle+1:h],:);
%      otherwise
%        emission = map;
%        transitions = dist;
%    end
%  end

  return;
end
