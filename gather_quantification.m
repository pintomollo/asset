function [datas, theta] = gather_quantification(mymovie, opts)

  npts = 500;

  nframes = size_data(mymovie.data);  
  theta = [-0.5:1/npts:0.5].';
  theta = theta(1:end-1);

  datas = NaN(nframes, npts);
  dist = NaN(nframes, 2);

  thresh = 1/10;

  type = 'direct';

  switch type
    case 'elliptic'
      theta = theta * 2*pi;
      range = [-pi pi];
    case 'linear'
      range = [-0.5 0.5];
    case 'direct'
      resolution = 0.5;
      tmp_pts = cell(nframes, 2);
  end

  for i=1:nframes
    
    warper = mymovie.data.warpers(i);

    if (isfield(mymovie.data.quantification, 'carth') & ~isempty(mymovie.data.quantification(i).carth))
      cortex = mymovie.data.quantification(i).carth;
    else
      cortex = mymovie.data.cortex(i).carth;
      if (opts.quantification.use_ruffles)
        cortex = insert_ruffles(cortex, mymovie.markers.ruffles(i).paths);
      end
    end


    %if (opts.quantification.use_ruffles)
    %  [ cortex, rescale] = insert_ruffles(cortex, mymovie.markers.ruffles(i).paths);
    %else
    %  rescale = false(size(cortex, 1), 1);
    %end

    cortex = carth2normalized(cortex, warper, opts);
    switch type
      case 'elliptic'
        pts = carth2elliptic(cortex, warper.reference.centers, warper.reference.axes_length, warper.reference.orientations);
        pts = pts(:, 1);
      case 'linear'
        %pts = carth2linear(cortex, mymovie.markers.ruffles(i).carth, mymovie.markers.ruffles(i).paths)
        pts = carth2linear(cortex);
      case 'direct'
        [pts, total_dist] = carth2linear(cortex);

        npts = length(pts);

        %ell_pts = carth2elliptic(cortex, warper.reference.centers, warper.reference.axes_length, warper.reference.orientations, 'radial');

        %thresh = median(diff(ell_pts(:,1))) * 1.5;

        goods = (cortex(:,1) > 0);
        max_pos = max(cortex(goods, 1));
        goods = (cortex(:,1) > (1-thresh)*max_pos);

        max_r = min(abs(cortex(goods,2)));
        post_indx = find(goods & abs(cortex(:,2)) == max_r, 1);

        goods = (cortex(:,1) < 0);
        max_pos = min(cortex(goods, 1));
        goods = (cortex(:,1) < (1-thresh)*max_pos);

        max_r = min(abs(cortex(goods,2)));
        ant_indx = find(goods & abs(cortex(:,2)) == max_r, 1);

        %goods = (abs(ell_pts(:,1)) <= thresh);
        %max_r = max(ell_pts(goods,2));

        %post_indx = find(goods & ell_pts(:,2) == max_r, 1);

        %goods = (abs(abs(ell_pts(:,1)) - pi) < thresh);
        %max_r = max(ell_pts(goods,2));
        %ant_indx = find(goods & ell_pts(:,2) == max_r, 1);

        pts = pts * total_dist;
        pts = pts - pts(post_indx);
        if (ant_indx > post_indx)
          pts = [pts(ant_indx:end, :) - total_dist; pts(1:ant_indx-1, :)];
        else
          pts = [pts(ant_indx:end, :); pts(1:ant_indx-1, :) + total_dist];
        end
        dist(i,1) = pts(1);
        dist(i,2) = total_dist + pts(1);

        tmp_pts{i,1} = pts;
        tmp_pts{i,2} = [ant_indx:npts 1:ant_indx-1];

        continue;
    end
    %full_theta = mymovie.data.warpers(i).warp(:,1);
    intens = mymovie.data.quantification(i).cortex;

    new_intens = interp_elliptic([pts, intens], theta, range);
    datas(i, :) = new_intens(:, 2);
    %bkg(i,:) = mymovie.data.quantification(i).bkg;
  end

  if (strncmp(type, 'direct', 6))
    %dist = dist * opts.pixel_size;
    boundaries = max(abs(dist));

    if (any(boundaries > 1e4))
      keyboard
      error('Too many elements');
    end

    boundaries(1) = -boundaries(1);
    full_indexes = [fliplr([0:-resolution:boundaries(1)]) resolution:resolution:boundaries(2)];
    npts = length(full_indexes);

    datas = NaN(nframes, npts);
    for i=1:nframes

      pts = tmp_pts{i,1}; % * opts.pixel_size;
      new_intens = interp_elliptic(pts, mymovie.data.quantification(i).cortex(tmp_pts{i,2}), full_indexes, dist(i,:));
      new_intens(full_indexes < dist(i,1) | full_indexes > dist(i,2),:) = NaN;
      datas(i, :) = new_intens(:, 2);
    end

    theta = full_indexes;
  end

  return;
end
