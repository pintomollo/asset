function [path, emission, transitions] = dynamic_programming(img, params, weight, weight_params, opts, retrieve_maps)

  if (nargout > 1 & (nargin < 6 | ~retrieve_maps))
    do_probs = true;
  else
    do_probs = false;
  end

  if (nargin < 6)
    retrieve_maps = false;
  end

  [h,w] = size(img);

  nhood = params.nhood;
  half = floor(nhood/2);

  switch params.dp_method
    case 'inverted'
      indxs = [h:-1:1];
    case 'backward'
      indxs = [1:h h-1:-1:1];
    case 'double'
      indxs = [1:h 1:h];
    otherwise
      indxs = 1:h;
  end
  nsteps = length(indxs);

  init = params.init; 
  path = NaN(h,1);

  dist = zeros(nsteps, w);
  map = zeros(nsteps, w);
  if (do_probs)
    emission = zeros(nsteps, nhood, w);
    transitions = zeros(1,nhood+2);
  end

  nstart = 0;
  is_circular = strncmp(params.prohibit, 'none', 4);
  spawn = params.spawn_percentile;

  if (isempty(weight_params))
    wimg = weight(img);
  else
    wimg = weight(img, weight_params);
  end

  no_spawn = all(~isfinite(img), 2);
  img(no_spawn, :) = mymean(img(:));

  if (any(isfinite(wimg(:))))
    wimg(no_spawn, :) = mymean(wimg(isfinite(wimg)));
  else
    wimg(no_spawn, :) = 1;
  end

  if (isempty(init))
    init = zeros(1,w);
  elseif (length(init) == 1)
    tmp = init;
    init = Inf(1,w);
    init(1,tmp) = 0;

    nstart = 1;
  else
    nstart = sum(~isinf(init));
  end

  spawn_gap = Inf;
  if (~isempty(spawn) & all(isfinite(spawn) & spawn <= 1 & spawn >= 0))
    switch params.spawn_type
      case 'full'
        spawn_gap = prctile(wimg(isfinite(wimg)), 100*spawn(1)) * (params.alpha) + (1-params.alpha)*spawn;
      case 'end'
        tmp_img = wimg(end-10:end, :);
        spawn_gap = prctile(tmp_img(isfinite(tmp_img)), 100*spawn(1)) * (params.alpha) + (1-params.alpha)*spawn(1);
    end
  end
  if (~isfinite(spawn_gap))
    spawn_gap = Inf;
    no_spawn(:) = true;
  end
  spawn_dist = init;

  prev_line = img(1,:);
  if (do_probs)
    [dist(1,:), map(1,:), emission(1,:,:), transitions(1,1)] = dp_score_mex([], prev_line, wimg(1,:), init, [], params, is_circular);
  else
    [dist(1,:), map(1,:)] = dp_score_mex([], prev_line, wimg(1,:), init, [], params, is_circular);
  end

  for j=2:nsteps

    indx = indxs(j);
    line = img(indx,:);
    
%    switch params.spawn_type
%      case 'current'
%        spawn_gap = prctile(wimg(indx, isfinite(wimg(indx, :))), 100*spawn) * (params.alpha) + (1-params.alpha)*spawn;
%      case 'previous'
%        tmp = wimg(1:indx, :);
%        spawn_gap = prctile(tmp(isfinite(tmp)), 100*spawn) * (params.alpha) + (1-params.alpha)*spawn;
%    end
    %if (~isfinite(spawn_gap))
    %  spawn_gap = 0;
    %end
    spawn_dist = spawn_dist + spawn_gap;

    if (do_probs)
      if (j==2)
        [dist(j,:), map(j,:), emission(j,:,:), transitions(1,2:end-1)] = dp_score_mex(prev_line, line, wimg(indx,:), dist(j-1,:), map(j-1,:), params, is_circular, spawn_dist);
      else
        [dist(j,:), map(j,:), emission(j,:,:)] = dp_score_mex(prev_line, line, wimg(indx,:), dist(j-1,:), map(j-1,:), params, is_circular, spawn_dist);
      end
    else
      if (no_spawn(indx))
        [dist(j,:), map(j,:)] = dp_score_mex(prev_line, line, wimg(indx,:), dist(j-1,:), map(j-1,:), params, is_circular);
      else
        [dist(j,:), map(j,:)] = dp_score_mex(prev_line, line, wimg(indx,:), dist(j-1,:), map(j-1,:), params, is_circular, spawn_dist);
      end
    end

    if (all(isinf(dist(j,:))))
      dist(j,:) = dist(j-1,:);
      map(j,:) = [1:w];
      line = prev_line;
    end

    prev_line = line;
  end

  if (do_probs)
    [junk, junk, junk, transitions(1,end)] = dp_score_mex([], prev_line, wimg(indxs(end),:), init, [], params, is_circular);
  end

  tmp = NaN(nsteps, 1);

  if (strncmp(params.spawn_type, 'best', 4))
    no_spawn = all(~isfinite(img), 2);

    if (spawn(1) < 1)
      spawn(1) = spawn(1) * w;
    end

    if (length(spawn) == 1)
      spawn = [spawn 0];
    else
      spawn(2) = spawn(2)*h;
    end

    tmp(end,1) = find(dist(end,:)==min(dist(end,:)),1);

    for j=1:nsteps-1
      if (map(end-j+1, tmp(end-j+1,1)) == 0)
        break;
      end

      tmp(end-j,1) = map(end-j+1,tmp(end-j+1,1));

      alt_start = find(dist(end-j,:)==min(dist(end-j,:)),1);
      if (abs(tmp(end-j, 1) - alt_start) > spawn(1) & h-j > spawn(2) & ~no_spawn(end-j))
        tmp(:,1) = NaN;
        tmp(end-j, 1) = alt_start;
      end
    end

  else

    if (isfield(params, 'final') & ~isempty(params.final))
      subdist = dist(end, params.final);
      tmp(end,1) = params.final(find(subdist == min(subdist), 1));
    else
      tmp(end,1) = find(dist(end,:)==min(dist(end,:)),1);
    end

    for j=1:nsteps-1
      if (map(end-j+1, tmp(end-j+1,1)) == 0)
        break;
      end

      tmp(end-j,1) = map(end-j+1,tmp(end-j+1,1));
    end
  end

  switch params.dp_method
    case 'inverted'
      path = tmp(indxs);
      if (do_probs)
        emission = emission(indxs, :, :);
      end
    case 'backward'
      path = tmp(end:-1:h,1);
      if (do_probs)
        emission = emission(end:-1:h, : , :);
      end
    case 'double'
      middle = floor(h/2);
      path = tmp([h+1:h+middle middle+1:h],1);

      if (do_probs)
        emission = emission([h+1:middle+h middle+1:h], :, :);
      end
    otherwise
      path = tmp;
  end

  if (params.force_circularity && abs(path(1,1)-path(end,1))>(nhood/2))
    
    check = [-half:half];

    if (nstart ~= 1)

      if (nstart < w & nstart > 0)
        dcheck = floor(nstart / 2) - 1;
        check = [-dcheck:dcheck];

        if (isempty(check))
          check = 0;
        end
      end
      check = check + path(end,1);
      check = check(check > 0 & check <= w);

      init = Inf(1,w);
      init(1, check) = 0;

      params.init = init;

      if (do_probs)
        [path, emission, transitions] = dynamic_programming(img, params, weight, weight_params, opts);
      else
        [path] = dynamic_programming(img, params, weight, weight_params, opts);
      end
    else
      check = check + path(1,1);
      check = check(check > 0 & check <= w);

      [junk, mins] = min(dist(1, check));

      path(end,1) = check(mins);

      for j=1:h-1
        new_path = map(end-j+1,path(end-j+1,1));
        if (path(end-j,1)==new_path)
          break;
        end
        path(end-j,1) = new_path;
      end
    end
  end

  if (retrieve_maps)
    switch params.dp_method
      case 'inverted'
        emission = map(indxs,:);
        transitions = dist(indxs,:);
      case 'backward'
        emission = map(end:-1:h,:);
        transitions = dist(end:-1:h,:);
      case 'double'
        middle = floor(h/2);
        emission = map([h+1:h+middle middle+1:h],:);
        transitions = dist([h+1:h+middle middle+1:h],:);
      otherwise
        emission = map;
        transitions = dist;
    end
  end

  return;
end
