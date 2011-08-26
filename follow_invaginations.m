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

  inner_thresh = 0.02;
  
  coef = -2;
  
  for nimg=1:nframes
    %nimg = floor(rand(1)*nframes)+1
    %nimg = 28;
    
    ruffles = mymovie.(type).ruffles(nimg).carth;
    nruffles = size(ruffles,1);
    paths = cell(nruffles, 1);

    if (nruffles > 0)

      inv_length = zeros(nruffles, 1);
      tmp_mins = zeros(nruffles, 1);
      bounds = mymovie.(type).ruffles(nimg).bounds;
      mins = bounds(:,1:2);
      maxs = bounds(:,3:4);

      img = double(load_data(mymovie.(channel), nimg));
      if (strncmp(type, 'markers',7))
        img = gaussian_mex(img, parameters.noise.gaussian);
        img = median_mex(img, parameters.noise.median);
        img = imnorm(img);
      end

      path = mymovie.(type).cortex(nimg).carth;

      polar_img = elliptic_coordinate(img, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg), parameters.safety);
      polar_size = size(polar_img);
      vert_indx = [1:polar_size(1)]';
      flipped_size = polar_size([2 1]);
      polar_img = polar_img.';

      dist_thresh = inner_thresh * polar_size(2);

      ell_path = carth2elliptic(path, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg));
      ell_path = elliptic2pixels(ell_path, polar_size, parameters.safety);
      ell_path = adapt_path(polar_size, ell_path);

      ell_ruffles = carth2elliptic(ruffles, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg));
      ell_ruffles = elliptic2pixels(ell_ruffles, polar_size, parameters.safety);

      ell_mins = carth2elliptic(mins, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg));
      ell_mins = elliptic2pixels(ell_mins, polar_size, parameters.safety);
      ell_maxs = carth2elliptic(maxs, centers(:,nimg), axes_length(:,nimg), orientations(1,nimg));
      ell_maxs = elliptic2pixels(ell_maxs, polar_size, parameters.safety);

      ell_mins = round(ell_mins(:, [2 1]));
      ell_mins(ell_mins(:,2) < 1, 2) = ell_mins(ell_mins(:,2) < 1, 2) + flipped_size(2);
      ell_mins(ell_mins(:,2) > flipped_size(2), 2) = ell_mins(ell_mins(:,2) > flipped_size(2), 2) - flipped_size(2);
      ell_maxs = round(ell_maxs(:, [2 1]));
      ell_maxs(ell_maxs(:,2) < 1, 2) = ell_maxs(ell_maxs(:,2) < 1, 2) + flipped_size(2);
      ell_maxs(ell_maxs(:,2) > flipped_size(2), 2) = ell_maxs(ell_maxs(:,2) > flipped_size(2), 2) - flipped_size(2);
      
      %ell_bounds = round([ell_mins ell_maxs]);
      %ell_bounds(ell_bounds < 1) = ell_bounds(ell_bounds < 1) + flipped_size(2);
      %ell_bounds(ell_bounds > flipped_size(2)) = ell_bounds(ell_bounds > flipped_size(2)) - flipped_size(2);

      [junk, map, dist] = dynamic_programming(polar_img, parameters.cortex_params, segment_func, parameters.cortex_weights, opts, true);

      dist_path = bilinear_mex(dist, vert_indx, ell_path);
      weight_img = segment_func(polar_img, parameters.cortex_weights);
      intens_path = bilinear_mex(weight_img, vert_indx, ell_path);
      
      mean_intens = mean(intens_path);
      thresh = mean_intens+ coef*std(intens_path);

      if (opts.verbosity == 3)
        figure;subplot(2,1,1);imshow(weight_img);hold on;
        plot([1:length(ell_path)],ell_path);
        subplot(2,1,2);hold on;plot([0 0; polar_size(2) polar_size(2)],[mean_intens thresh; mean_intens thresh])
      end

      for i=1:nruffles
        if (ell_mins(i,2) > ell_maxs(i,2))
          circular = true;
          portion = dist_path([ell_mins(i,2):end 1:ell_maxs(i,2)]);
        else
          circular = false;
          portion = dist_path([ell_mins(i,2):ell_maxs(i,2)]);
        end
        start = find(min(portion) == portion) + ell_mins(i,2) - 1;
        if (circular)
          start(start > flipped_size(2)) = start(start > flipped_size(2)) - flipped_size(2);
        end

        if (length(start) > 1)
          inv_dist = abs(start - ell_ruffles(i,1));
          start = start(find(inv_dist == min(inv_dist),1));
        end

        inv_path = backtrack(map, round(ell_path(start)), start);

       %   if (nimg == 120)
       %     myplot(inv_path)
       %     keyboard
       %   end

        inv_intens = weight_img(sub2ind(flipped_size, inv_path(:,1), inv_path(:,2)));
        worst = find([inv_intens(1:end-1); Inf] >= mean_intens, 1);
        indx = find(inv_intens(1:worst) < thresh, 1, 'last');
        
        inv_dist = ell_path(inv_path(:,2)) -inv_path(:,1);
        %min_indx = find(inv_dist == min(inv_dist(inv_dist > 0)), 1, 'last');
        last_neg = find(inv_dist <= 0, 1, 'last');
        if (isempty(last_neg))
          last_neg = 1;
        end

        if (last_neg <= indx)
          min_indx = find([true; (inv_dist(last_neg+1:indx) < inv_dist(last_neg:indx-1))], 1, 'last') + last_neg-1;

          %keyboard
          if (opts.verbosity == 3)
            subplot(2,1,2);
            plot(inv_path(:,1), inv_intens);
            subplot(2,1,1);
            myplot(inv_path(:,[2 1]),'g');
          end

          if (indx > 1 & inv_dist(indx) > dist_thresh)

            if (min_indx ~= 1)
              min_indx = min_indx - 1;
              inv_path(min_indx,1) = ell_path(inv_path(min_indx,2));
            end

            inv_path = inv_path(min_indx:indx,[2 1]);
            
            if (opts.verbosity == 3)
              myplot(inv_path,'r');
            end

            %keyboard
            
            inv_path = pixels2elliptic(inv_path,polar_size,parameters.safety);
            inv_path = elliptic2carth(inv_path,centers(:,nimg),axes_length(:,nimg),orientations(1,nimg));
            paths{i} = inv_path;
            inv_length(i,1) = sum(sqrt(sum(diff(inv_path, [], 1).^2, 2)));

            tmp_mins(i) = min_indx;

            is_valid = true;
            for j=1:i-1
              if (is_valid & inv_length(j,1) > 0 & all(paths{i}(end,:) == paths{j}(end,:)))
                if (tmp_mins(i) < tmp_mins(j))
                  paths{j} = [];
                  inv_length(j,1) = 0;
                else
                  paths{i} = [];
                  inv_length(i,1) = 0;
                  is_valid = false;
                end
              end
            end
          end
        end
      end
      mymovie.(type).ruffles(nimg).properties = [mymovie.(type).ruffles(nimg).properties inv_length];

      %keyboard
    end

    %keyboard
    mymovie.(type).ruffles(nimg).paths = paths;
    %keyboard
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
