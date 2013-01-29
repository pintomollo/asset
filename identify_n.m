function nindep = identify_n(file)

  if (ischar(file))
    data = load(file);
  else
    data = file;
  end

  if (isstruct(data))
    if (isfield(data, 'mymovie'))
      [fraction, max_width, cell_width, domain, pos] = domain_expansion(data.mymovie, data.opts);
      times = get_manual_timing(data.mymovie, data.opts);
      ground_truth = domain(1:times(end), :);
    else
      ground_truth = data.ground_truth;
    end
  else
    ground_truth = data;
  end
  [sh,sw,sn] = size(ground_truth);

  nindep = NaN(sn, 4);

  for p = 1:sn
    pacf = NaN(sh+sw, 21);
    bound = NaN(sh+sw, 2);
    for h = 1:sh
      signal = ground_truth(h,:,p);
      goods = isfinite(signal);
      if (any(goods))
        [pacf(h,:), junk, bound(h,:)] = parcorr(signal(goods));
      end
    end
    for w = 1:sw
      signal = ground_truth(:,w,p);
      goods = isfinite(signal);
      if (any(goods))
        [pacf(sh+w,:), junk, bound(sh+w,:)] = parcorr(signal(goods));
      end
      %[pacfw(1:ngoods+1, w), junk, bounds2(h,:)] = parcorr(signal(isfinite(signal)));
    end

    rescaled = bsxfun(@rdivide, pacf, bound(:,1));

    means = mymean(rescaled(1:sh,:), 1);
    nindep(p,1) = find(means>=1, 1, 'last') + 1;

    means = mymean(pacf(1:sh,:));
    nindep(p,2) = find(abs(means)>=mymean(bound(1:sh,1)), 1, 'last') + 1;

    means = mymean(rescaled(sh+1:end,:));
    nindep(p,3) = find(means>=1, 1, 'last') + 1;

    means = mymean(pacf(sh+1:end,:));
    nindep(p,4) = find(abs(means)>=mymean(bound(sh+1:end,1)), 1, 'last') + 1;
  end


  return;

  keyboard

  means = mean(rescaled(sh+1:end,:));

  bounds = [mean(bound(1:sh,:)); mean(bound(sh+1:end,:))];

  means = [means; mean(bsxfun(@ge, pacf, bound(:,1)) | bsxfun(@le, pacf, bound(:,2)))];
  bounds = [bounds; [0.5 -0.5]];

  means = [means; mean(abs(pacf))];
  bounds = [bounds; mean(bound)];

  means = [means; mean(bsxfun(@rdivide, pacf, bound(:,1)))];
  bounds = [bounds; [1 -1]];

  keyboard


  orig = inpaint_nans(ground_truth);

  min_width = 21;
  nsteps = 20;

  widths = linspace(min_width, sw, nsteps);
  widths = round(widths);
  widths = widths + (1 - mod(widths, 2));
  end_pos = pos(end);

  for w=widths
    ground_truth = imresize(orig, [sh w]);
    pos = [1:w] - ((w-1)/2) - 1;
    pos = pos * (end_pos / pos(end));

    name = [mymovie.experiment num2str(w) '.mat'];
    save(name, 'ground_truth', 'pos');

    find_kymograph(name, 'fit_full', true,  'config_fitting', 'param_distribution', 'integrate_sigma', true, 'max_iter', 2000);
  end
  
%{
  log2(size(domain))
  for i=1:nlogs
    img = imresize(img);
    save()
    find_kymograph()
  end
%}
  return;
end
