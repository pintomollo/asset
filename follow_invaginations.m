function mymovie = follow_invaginations(mymovie, opts)

  if (strncmp(opts.segmentation_type, 'markers', 7) & isfield(mymovie, 'markers') & ~isempty(mymovie.markers))
    type = 'markers';
    channel = 'cortex';
    [nframes, imgsize] = size_data(mymovie.cortex);
    segment_func = @intens;
  else
    type = 'dic';
    channel = 'dic';
    [nframes, imgsize] = size_data(mymovie.dic);
    segment_func = @edges_intens;
  end

  parameters = opts.segmentation_parameters.(type);
  centers = mymovie.(type).centers;
  axes_length = mymovie.(type).axes_length;
  orientations = mymovie.(type).orientations;
  neighbors = mymovie.(type).neighbors;

  %inner_thresh = 0.02;
  jump_coef = 0.05;
  coef = 3.5;
  straighten_coef = 0.5;
  bound_thresh = 50;

  for i = 1:nframes
    nimg = i;
    %nimg = floor(rand(1)*nframes)+1
    %nimg = 16 + i

    ruffles = mymovie.(type).ruffles(nimg).carth;
    valid_ruffles = all(~isnan(ruffles), 2);
    nruffles = size(ruffles,1);

    paths = cell(nruffles, 1);
    inv_length = NaN(nruffles, 1);
    mins = NaN(nruffles, 2);
    maxs = mins;

    if (nruffles > 0)

      bounds = mymovie.(type).ruffles(nimg).bounds;
      nbounds = size(bounds, 1);
      mins(1:nbounds, :) = bounds(:,1:2);
      maxs(1:nbounds, :) = bounds(:,3:4);

      img = double(load_data(mymovie.(channel), nimg));
      img = mask_neighbors(img, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), neighbors(nimg), opts);

      if (strncmp(type, 'markers',7))
        img = gaussian_mex(img, parameters.noise.gaussian);
        img = median_mex(img, parameters.noise.median);
        img = imnorm(img);
      end

      path = mymovie.(type).cortex(nimg).carth;

      polar_img = elliptic_coordinate(img, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), parameters.safety);
      polar_size = size(polar_img);
      polar_size = polar_size(:);
      vert_indx = [1:polar_size(1)]';
      flipped_size = polar_size([2 1]);
      polar_img = polar_img.';

      %dist_thresh = inner_thresh * polar_size(2);
      jump_thresh = jump_coef * polar_size(2);

      ell_path = carth2elliptic(path, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg));
      ell_path = elliptic2pixels(ell_path, polar_size, axes_length(:, nimg), parameters.safety);
      ell_path = adapt_path(polar_size, ell_path);

      ell_ruffles = carth2elliptic(ruffles, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg));
      ell_ruffles = elliptic2pixels(ell_ruffles, polar_size, axes_length(:, nimg), parameters.safety);

      ell_mins = carth2elliptic(mins, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg));
      if (isempty(ell_mins))
        ell_mins = mins;
      else
        ell_mins = elliptic2pixels(ell_mins, polar_size, axes_length(:, nimg), parameters.safety);
        ell_mins = round(ell_mins(:, [2 1]));
        ell_mins(ell_mins(:,2) < 1, 2) = ell_mins(ell_mins(:,2) < 1, 2) + flipped_size(2);
        ell_mins(ell_mins(:,2) > flipped_size(2), 2) = ell_mins(ell_mins(:,2) > flipped_size(2), 2) - flipped_size(2);
      end

      ell_maxs = carth2elliptic(maxs, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg));
      if (isempty(ell_maxs))
        ell_maxs = maxs;
      else
        ell_maxs = elliptic2pixels(ell_maxs, polar_size, axes_length(:, nimg), parameters.safety);
        ell_maxs = round(ell_maxs(:, [2 1]));
        ell_maxs(ell_maxs(:,2) < 1, 2) = ell_maxs(ell_maxs(:,2) < 1, 2) + flipped_size(2);
        ell_maxs(ell_maxs(:,2) > flipped_size(2), 2) = ell_maxs(ell_maxs(:,2) > flipped_size(2), 2) - flipped_size(2);
      end

      %ell_bounds = round([ell_mins ell_maxs]);
      %ell_bounds(ell_bounds < 1) = ell_bounds(ell_bounds < 1) + flipped_size(2);
      %ell_bounds(ell_bounds > flipped_size(2)) = ell_bounds(ell_bounds > flipped_size(2)) - flipped_size(2);

      parameters.cortex_params.beta = (1 - straighten_coef)*parameters.cortex_params.beta + straighten_coef;
      [junk, map, dist] = dynamic_programming(polar_img, parameters.cortex_params, segment_func, parameters.cortex_weights, opts, true);

      dist_path = bilinear_mex(dist, vert_indx, ell_path);
      weight_img = segment_func(polar_img, parameters.cortex_weights);
      %intens_path = bilinear_mex(weight_img, vert_indx, ell_path);
      %intens_bkg = bilinear_mex(weight_img, vert_indx, ell_path/2);
      %intens_path = bilinear_mex(polar_img, vert_indx, ell_path);

      %mean_intens = prctile(intens_path, 50);
      %mean_bkg = prctile(intens_bkg, 50);
      %thresh = mean_intens+ coef*std(intens_path);
      %thresh = mean_intens + coef*(mean_bkg - mean_intens);
      %thresh = prctile(weight_img, coef, 2);
      thresh = prctile(weight_img, 50, 2) - coef*std(weight_img, [], 2);

      if (opts.verbosity == 3)
        figure;subplot(2,1,1);imshow(weight_img);hold on;
        plot([1:length(ell_path)],ell_path);
        %subplot(2,1,2);hold on;plot([0 0; polar_size(2) polar_size(2)],[mean_intens thresh; mean_intens thresh])
        subplot(2,1,2);hold on;plot(thresh, 'r');
      end

      for r=1:nruffles
        if (~valid_ruffles(r))
          continue;
        end

        init_pos = ell_mins(r,2);
        if (ell_mins(r,2) < ell_maxs(r,2))
          circular = false;
          portion = dist_path([ell_mins(r,2):ell_maxs(r,2)]);
        elseif (isnan(init_pos))
          bound_size = (polar_size(2)/bound_thresh);
          curr_bounds = round([ell_ruffles(r, 2) - bound_size ell_ruffles(r, 2) + bound_size]);
          if (curr_bounds(1) < 1)
            curr_bounds(1) = curr_bounds(1) + polar_size(1);
          end
          if (curr_bounds(2) > polar_size(1))
            curr_bounds(2) = curr_bounds(2) - polar_size(1);
          end

          if (curr_bounds(1) < curr_bounds(2))
            circular = false;
            portion = dist_path([curr_bounds(1):curr_bounds(2)]);
          else
            circular = true;
            portion = dist_path([curr_bounds(1):end 1:curr_bounds(2)]);
          end

          init_pos = curr_bounds(1);
        else
          circular = true;
          portion = dist_path([ell_mins(r,2):end 1:ell_maxs(r,2)]);
        end
        start = find(min(portion) == portion) + init_pos - 1;
        if (circular)
          start(start > flipped_size(2)) = start(start > flipped_size(2)) - flipped_size(2);
        end

        if (length(start) > 1)
          inv_dist = abs(start - ell_ruffles(r,1));
          start = start(find(inv_dist == min(inv_dist),1));
        elseif (length(start) == 0)
          inv_length(r,1) = 0;
          continue;
        end

        inv_path = backtrack(map, round(ell_path(start)), start);
        inv_intens = weight_img(sub2ind(flipped_size, inv_path(:,1), inv_path(:,2)));

        val_test = ([inv_intens(1:end-1); Inf] >= thresh(inv_path(1,1):-1:1));
        %transitions = find([2; diff(val_test)]);
        %domains = diff([transitions; size(inv_intens, 1)]);
        %domain_vals = val_test(transitions);
        [transitions, domains, domain_vals] = boolean_domains(val_test);

        good_domains = ((~domain_vals & (domains > [0; domains(1:end-1)])) | (domain_vals & domains < jump_thresh));
        indx = find(~good_domains, 1);

        if (isempty(indx))
          worst = length(inv_intens);

          %%%%%%%%%%%%% Reached the end, should follow them at pos + size/2 backward
          %%%%%%%%%%%%% so re-run the dyn_prog backward, store it.
          %%%%%%%%%%%%% Need also to check if another path reaches there
        else
          worst = transitions(indx);
        end

        if (opts.verbosity == 3)
          subplot(2,1,2);
          plot(inv_path(:,1), inv_intens);
          subplot(2,1,1);
          myplot(inv_path(:,[2 1]),'g');
        end

        if (worst > 1)
          indx = find(inv_intens(1:worst) < thresh(worst:-1:1), 1, 'last');
          if (isempty(indx))
            indx = 1;
          end

          inv_dist = ell_path(inv_path(:,2)) -inv_path(:,1);
          last_neg = find(inv_dist <= 0, 1, 'last');
          if (isempty(last_neg))
            last_neg = 1;
          end

          %if (last_neg+1 < indx & inv_dist(indx) > dist_thresh)
          if (last_neg+1 < indx)

            inv_path = inv_path(last_neg+1:indx,[2 1]);
            inv_path = pixels2elliptic(inv_path,polar_size, axes_length(:,nimg),parameters.safety);
            inv_path = elliptic2carth(inv_path,centers(:,nimg),axes_length(:,nimg),orientations(1,nimg));
            paths{r} = inv_path;
            inv_length(r,1) = sum(sqrt(sum(diff(inv_path, [], 1).^2, 2)));

            is_valid = true;
            for j=1:r-1
              if (is_valid && inv_length(j,1) > 0 && all(paths{r}(end,:) == paths{j}(end,:)))
                if (sum((paths{r}(1,:) - ruffles(r,:)).^2) < sum((paths{j}(1,:) - ruffles(j,:)).^2))
                  paths{j} = [];
                  inv_length(j,1) = 0;
                else
                  paths{r} = [];
                  inv_length(r,1) = 0;
                  is_valid = false;
                end
              end
            end
          end
        end
      end
      mymovie.(type).ruffles(nimg).properties = [mymovie.(type).ruffles(nimg).properties inv_length];

      if (opts.verbosity == 3)
        subplot(2,1,1);
        for p=1:length(paths)
          if (~isempty(paths{p}))
            inv_path = carth2elliptic(paths{p},centers(:,nimg),axes_length(:,nimg),orientations(1,nimg));
            inv_path = elliptic2pixels(inv_path,polar_size, axes_length(:,nimg),parameters.safety);
            myplot(inv_path,'r');
          end
        end
      end
    end

    mymovie.(type).ruffles(nimg).paths = paths;
  end

  return;
end

function path = backtrack(map, row, col)

  path = zeros(row, 2);

  path(1,:) = [row, col];
  for i=row-1:-1:1
    path(end-i+1,:) = [i, map(i+1, path(end-i,2))];
  end

  return;
end

function path = forwardtrack(map, dist, col)

  path = zeros(size(map, 1), 1);

  path(1) = col;
  for i=2:length(path)
    indexes = find(map(i, :) == path(i-1));
    best = find(dist(i, indexes) == min(dist(i, indexes)));
    path(i) =  indexes(best);
  end

  return;
end
