function [p] = fit_split(mymovie, opts)

  %P = [  0.3676 0.049336 0.011196 0.70123 0.1047]
  %P = [0.3676 0.2 0.011196 0.70123 0.1047]
  %P = [0.3676 0.05 0.05 0.45 0.15]
  %P = [0.9809 0.0494 0.1849 0.4497 0.1573]

  %P(end) = P(end)*100;

  [nframes, ssize] = size_data(mymovie.dic);
  estim_only = false; 
  ellipses = cell(nframes, 1);
 
  if (false)
  for i=1:nframes
    img = imnorm(double(load_data(mymovie.dic, i)));

    imshow(img);hold on;

    img = imadm_mex(img);
    thresh = graythresh(img);
    img = (img > thresh*0.5*(max(img(:))) );

    [ellipse, estim] =  split_cells(img, estim_only, opts, P);

        for j = 1:size(ellipse, 1)
          draw_ellipse(ellipse(j, 1:2), ellipse(j, 3:4), ellipse(j, 5));
        end
        hold off;

    ellipses{i} = ellipse;

    print('-dpng', ['./PNG/separate1-' mymovie.experiment  '-' num2str(i) '.png']);
  end

  save([mymovie.experiment '-split.mat'], 'ellipses');
  else
  load([mymovie.experiment '-split.mat']);
  end
 
  [n,m] = size(ellipses{1});
  real_ell = NaN(n, m, nframes);
  real_ell(:,:,1) = ellipses{1};

  dist_thresh = 50^2;

  for i=2:nframes
    tmp_ell = ellipses{i};
    for k=1:size(tmp_ell, 1)
      found = false;
      mpos = mymean(real_ell(:,1:2,:), 3);
      for j=1:size(mpos, 1)
        if (sum((mpos(j, :) - tmp_ell(k, 1:2)).^2) < dist_thresh)
          found = true;
          real_ell(j, :, i) = tmp_ell(k, :);
          break;
        end
      end

      if (~found)
        real_ell(end+1,:,:) = NaN;
        real_ell(end, :, 1) = tmp_ell(k,:);
      end
    end
  end

  goods = sum(~isnan(real_ell),3);
  goods = (goods(:,1) > nframes/2);

  %avg_ell = mymean(real_ell(goods, :, :), 3);
  avg_ell = NaN(0,m);

  img = NaN(ssize);
  img(1:end) = imnorm(double(load_data(mymovie.dic, 1)));

  figure;
  imshow(imnorm(double(load_data(mymovie.dic, 1))));
  hold on;
  for i=1:size(real_ell, 1)
    if (goods(i))
      tmp = real_ell(i, :, :);
      tmp = reshape(tmp(~isnan(tmp)), 5, []);
      avg_ell(end+1,:) = median(tmp, 2);

      draw_ellipse(avg_ell(end,:));
    end
  end
  drawnow
  avg_area = prod(avg_ell(:,3:4), 2)*pi;

  opt = cmaes('defaults');
  opt.MaxFunEvals = 10000;
  opt.TolFun = 1e-5;
  opt.SaveFilename = '';
  opt.LogFilenamePrefix = ['fit-split' mymovie.experiment '_'];
  opt.EvalParallel = 'yes';
  opt.LogPlot = 0;
  opt.Restarts = 1;

  p0 = [0.5 0.2 0.05 0.4 0.15].';
  nparams = length(p0);
  
  display(['IC (' num2str(p0.') ')']);
  [p, fval, ncoutns, stopflag, out] = cmaes(@error_function, p0, 0.1, opt);
  display(['Best (' num2str(p.') ')']);

  return;

  function err_all = error_function(p_all)

    [curr_nparams, nevals] = size(p_all);
    flip = false;
    if (nevals == nparams)
      flip = true;
      p_all = p_all.';
      nevals = curr_nparams; 
    end
    err_all = zeros(1, nevals);


    for nimg=1:nframes
      img(1:end) = imnorm(double(load_data(mymovie.dic, nimg)));

      dic_img = img;

      img(1:end) = imadm_mex(img);
      thresh = graythresh(img);
      img(1:end) = (img > thresh*0.5*(max(img(:))) );
  
      for i = 1:nevals
        match_ell = false(size(avg_ell, 1), 1);

        new_p = p_all(:, i).';
        new_p(end) = new_p(end)*100;

        [ellipse] = split_cells(img, estim_only, opts, new_p);

        %if (nevals > 1)
        %figure;imshow(dic_img);hold on;

        %  for j=1:length(match_ell)
        %  draw_ellipse(avg_ell(j, 1:2), avg_ell(j, 3:4), avg_ell(j, 5), 'r');
        %  end
        %end

        for k=1:size(ellipse, 1)
          %if (nevals > 1)
          %draw_ellipse(ellipse(k, 1:2), ellipse(k, 3:4), ellipse(k, 5));
          %end

          ell_area = prod(ellipse(k,3:4))*pi;
          max_overlap = Inf;
          for j=1:size(avg_ell, 1)
            curr_overlap = ellipse_ellipse_area_mex(ellipse(k,5), ellipse(k,3), ellipse(k,4), ellipse(k,1), ellipse(k,2), avg_ell(j,5), avg_ell(j,3), avg_ell(j,4), avg_ell(j,1), avg_ell(j,2));

            if (sum((avg_ell(j, 1:2) - ellipse(k, 1:2)).^2) < dist_thresh)
              match_ell(j) = true;
              max_overlap = avg_area(j) + ell_area - 2*curr_overlap;

              break;
            elseif (curr_overlap < max_overlap)
              max_overlap = (ell_area - curr_overlap);
            end
          end

          err_all(i) = err_all(i) + max_overlap;
        end

        err_all(i) = err_all(i) + 1.5*sum(avg_area(~match_ell));
      end

      %if (nevals>1)
      %keyboard
      %end
    end

    for i = 1:nevals
      new_p = p_all(:, i).';

      bads = (new_p < 0);
      if (any(bads))
        err_all(i) = err_all(i) + sum(exp(-10*new_p(bads)));
      end
    end

    if (flip)
      err_all = err_all.';
    end

    if (all(isinf(err_all)))
      % ML crashes when all values are non-numerical
      err_all(1) = 1e5;
    end

    return;
  end
end
