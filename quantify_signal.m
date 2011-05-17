function signal = quantify_signal(imgs, positions, opts)

  if (nargin == 2)
    if (isstruct(positions))
      opts = positions;
      positions = {};
    else
      opts = get_struct('ASSET');
    end
  elseif (nargin == 1)
    opts = get_struct('ASSET');
    positions = {};
  end

  channel_name = opts.quantification.channel;
  field_name = opts.quantification.field;
  safety = 1.2;

  if (isstruct(imgs))
    mymovie = imgs;
    imgs = [];

    if (opts.recompute | ~isfield(mymovie.(channel_name), 'eggshell') | isempty(mymovie.(channel_name).eggshell))
    
      mymovie = duplicate_segmentation(mymovie, channel_name, opts);
    end
    
    nframes = size_data(mymovie.(channel_name));
  else
    mymovie = [];
    nframes = size(imgs, 3);
  end

  window_shape = opts.quantification.window_shape;

  if (ischar(window_shape))
    window_size = opts.quantification.window_size / opts.pixel_size;

    switch window_shape
      case 'gaussian'
        window_shape = fspecial(window_shape, ceil(window_size), opts.quantification.window_params / opts.pixel_size);
        
      otherwise
        window_shape = fspecial(window_shape, ceil(window_size));
    end

    window_shape = window_shape / sum(window_shape(:));
  end

  %selem = strel('disk', ceil(opts.quantification.norm_shift / opts.pixel_size));
  pixel_shift = (opts.quantification.norm_shift / opts.pixel_size);

  for i=1:nframes

    %nimg = i
    nimg = randi(nframes)
    %nimg = 62;

    if (isempty(mymovie))
      img = imgs(:,:,nimg);
      pos = positions{nimg};
    else
      img = double(load_data(mymovie.(channel_name), nimg));
      pos = mymovie.(channel_name).(field_name)(nimg).carth;
    end

    img = imfilter(img, window_shape);
    %values = bilinear(img, pos(:,1), pos(:,2));

    switch opts.quantification.normalize
      case 'cytoplasm'
        [center, axes_length, orientation] = fit_ellipse(pos);

        polar_img = elliptic_coordinate(img, center, axes_length, orientation, safety);
        polar_size = size(polar_img);
        %polar_size = polar_size([2 1]);
        egg_path = carth2elliptic(pos, center, axes_length, orientation);
        
        correction = sqrt(cos(egg_path(:,1)).^2 + (axes_length(2)/axes_length(1) * sin(egg_path(:,1))).^2);

        egg_path = elliptic2pixels(egg_path, polar_size, safety);

        inside = [egg_path(:,1) egg_path(:,2) - correction * pixel_shift];
        %inside = elliptic2pixels(inside, polar_size, safety);

        values = bilinear(polar_img, egg_path(:,2), egg_path(:,1), true);
        cyto = bilinear(polar_img, inside(:,2), inside(:,1), true);
        %mask = roipoly(size(img,1),size(img,2), pos(:,1), pos(:,2));
        %mask = imerode(mask, selem);
        %inside = bwboundaries(mask);
        %inside = inside{1};
        %ncyto = size(inside, 1);
        %npts = size(pos, 1);
        %cyto = bilinear(img, inside(:,2), inside(:,1));
        %cyto = interp1q([1:ncyto].', cyto, [1:(ncyto-1)/(npts-1):ncyto].');

        %if (i > 50)
       % 
       % figure;imshow(imnorm(polar_img));
       % hold on;
       % myplot(egg_path(:,[2 1]));
       % myplot(inside(:,[2 1]), 'r');

        %figure;
        %plot(values);
        %hold on;
        %plot(cyto, 'r');

        %keyboard

        %values = values - cyto; 
        %end
      %case 'background'
      case 'perpendicular'
        
        ell_cortex = abs(carth2elliptic(pos, mymovie.data.centers(:, nimg), mymovie.data.axes_length(:, nimg), mymovie.data.orientations(1, nimg)));
        indx = find(ell_cortex(:,1)==min(ell_cortex(:,1)), 1);

        pos = pos([indx:end 1:indx-1], :);
        cortex = pos;        

        if (~isempty(mymovie))
          mymovie.(channel_name).(field_name)(nimg).carth = pos;
        end

        derivative_order = 2;

        % smoothing the "signal"
        switch derivative_order
          case 2
            cortex = [cortex(end,:); cortex; cortex(1,:)];
            dpos = (cortex(3:end, :) - cortex(1:end-2, :)) / 2;
          case 3
            cortex = [cortex(end-1:end, :); cortex; cortex(1:2, :)];
            dpos = (cortex(1:end-4, :) - 8*cortex(2:end-3, :) + 8*cortex(4:end-1, :) - cortex(5:end, :)) / 12;
          otherwise
            cortex = [cortex; cortex(1,:)];
            dpos = diff(cortex, [], 1);
        end
        dperp = [-dpos(:,2) dpos(:,1)];
        dperp = bsxfun(@rdivide, dperp, hypot(dperp(:,1), dperp(:,2))); 

        nbins = 16;

        dpos = [-nbins:nbins]*2;
        nbins = length(dpos);

        all_pos_x = bsxfun(@plus, dperp(:,1) * dpos, pos(:,1));
        all_pos_y = bsxfun(@plus, dperp(:,2) * dpos, pos(:,2));

        all_pos_x = all_pos_x(:);
        all_pos_y = all_pos_y(:);
        values = bilinear(img, all_pos_x, all_pos_y);

        values = reshape(values, [size(pos,1) nbins]);

        %figure;imagesc(values)

        all_pos_x = reshape(all_pos_x, [size(pos,1) nbins]); 
        all_pos_y = reshape(all_pos_y, [size(pos,1) nbins]); 

        figure;imshow(imnorm(img));
        hold on;
        %quiver(pos(:,1), pos(:, 2), dperp(:,1), dperp(:,2), 'r');
        %plot(all_pos_x(:,[1 end]), all_pos_y(:,[1 end]), 'g')


        %figure;imagesc(values);
        params = fit_signal(values, dpos);
        mymovie.(channel_name).quantification(nimg).front = params;

        keyboard

        new_values = extract_ridge(params, pos, dperp, opts);
        %figure;imagesc(front_function(values, dpos));
        %keyboard

        values = new_values;

      otherwise
        values = bilinear(img, pos(:,1), pos(:,2));
    end
    
    if (isempty(imgs))
      mymovie.(channel_name).quantification(nimg).(field_name) = values;
    else
      signal{nimg} = values;
    end
  end

  if (isempty(imgs))
    signal = mymovie;
  end

  return;
end

function best = fit_signal(values, pos)

  max_val = max(values(:));

  wsize = floor(length(pos) / 3);
  best = NaN(size(values, 1), 6);
  lbound = [-wsize/4 1 0 0 -0.1 0];
  ubound = [wsize/4 wsize max_val max_val 0.1 max_val];

  range = ubound - lbound;
  factor = 100;
    
  % Even that few iterations seems to be alright
  %opts = optimset('Display','off', 'Algorithm', 'levenberg-marquardt', 'MaxFunEval', 100);
  opts = optimset('Display','off', 'Algorithm', 'levenberg-marquardt');

  for i=1:size(values, 1)
    signal = values(i, :);
    estim = estimate_front(signal, pos);

    %if (any(isinf(estim) | isnan(estim)))
    %  keyboard;
    %end
    
    best(i, :) = lsqnonlin(@fit_front, estim, [], [], opts);
    %best(i, :) = lsqcurvefit(@fit_front, estim, pos, signal, lbound, [], opts);
  end

  return;

  function err = fit_front(params)

    estim_front = (atan((pos - params(1)) * params(5))/pi + 0.5) * params(6) + ...
                  exp(-((pos - params(1)).^2) / (2*(params(2)^2)))*params(3) + params(4);
    
    %estim_front = front_function(params, pos);

    err = abs(estim_front - signal);
    penalty = exp(factor*sum((lbound(params < lbound) - params(params < lbound)) ./ range(params < lbound))) + exp(factor*sum((params(params > ubound) - ubound(params > ubound)) ./ range(params > ubound))) + exp(factor*sum(abs(imag(params)))) - 3;

    if (isnan(penalty) | isinf(penalty))
      penalty = 1e10;
    end
    err = err + penalty;
  
    return;
  end
end
