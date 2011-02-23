function mymovie = detect_cytokinesis(mymovie, opts)

  %keyboard

  min_thresh = 10;
  thresh = 5;
  angl_thresh = 0.2;

  if (strncmp(opts.segmentation_type, 'markers', 7) & isfield(mymovie, 'markers') & ~isempty(mymovie.markers))
    type = 'markers';
    [nframes, imgsize] = size_data(mymovie.cortex);
  else
    type = 'dic';
    [nframes, imgsize] = size_data(mymovie.dic);
  end

  if (~isfield(mymovie.(type), 'ruffles'))
    mymovie = track_ruffles(mymovie, opts);
  end

  if (isfield(mymovie.(type), 'inverted'))
    actu_invert = mymovie.(type).inverted;
  else
    actu_invert = false;
  end
  
  first = true;
  cyto_ok = false;

  cyt1 = -1;
  cyt2 = -1;

  %keyboard

  for i=nframes:-1:1
    
    ell_pts = carth2elliptic(mymovie.(type).cortex(i).carth, mymovie.(type).centers(:,i),mymovie.(type).axes_length(:,i),mymovie.(type).orientations(1,i));
    %myplot(ell_pts);


    depth = mymovie.(type).ruffles(i).properties(:,1);
    %a(i) = max(depth);

    if (first)
      %max_indx = find(depth == max(depth));

      ell_pos = carth2elliptic(mymovie.(type).ruffles(i).carth, mymovie.(type).centers(:,i),mymovie.(type).axes_length(:,i),mymovie.(type).orientations(1,i));

      if (isempty(ell_pos))
        continue;
      end

      %keyboard

      dists = abs(ell_pos(:,1) - pi/2 - pi*(ell_pos(:,1) > pi));
      max_indx = find(dists == min(dists) & dists < angl_thresh);
      if (isempty(max_indx))
          continue;
      end
      
      angl = ell_pos(max_indx, 1); 

      if ((angl < pi & angl > pi/2) | (angl > pi & angl < 3*pi/2))
        disp(['A-P were inverted !']);
        inverted = ~actu_invert;
      else
        inverted = actu_invert;
      end

      fields = fieldnames(mymovie);
      for f = 1:length(fields)
        if (~isempty(mymovie.(fields{f})) && isfield(mymovie.(fields{f}), 'orientations') && ~isempty(mymovie.(fields{f}).orientations))
          if ((~isfield(mymovie.(fields{f}), 'inverted') & inverted) | (isfield(mymovie.(fields{f}), 'inverted') & mymovie.(fields{f}).inverted ~= inverted))
            tmp_orient = mymovie.(fields{f}).orientations;
            tmp_orient = tmp_orient + pi;
            tmp_orient = tmp_orient - 2*pi*(tmp_orient > 2*pi);
            mymovie.(fields{f}).orientations = tmp_orient;
            mymovie.(fields{f}).inverted = inverted;
          end
        end
      end

      first = false;

      cyt1 = mymovie.(type).ruffles(i).cluster(max_indx,1);
      %values(:,i) = junk(end-1:end).';
      if (length(depth) == 1)
        continue;
      end

      dist = abs(ell_pos(:,1)-(2*pi - angl));
      second_indx = find(dist == min(dist));
      cyt2 = mymovie.(type).ruffles(i).cluster(second_indx,1);
    else
      if (cyt2 > 0)
        cyt2 = mymovie.(type).ruffles(i).cluster(cyt2,1);
      elseif (cyt2 == -1)
        if (cyt1 == 0)
          cyt2 = 0;
        elseif (size(mymovie.(type).ruffles(i).properties,1) > 1)
          ell_pos = carth2elliptic(mymovie.(type).ruffles(i).carth, mymovie.(type).centers(:,i),mymovie.(type).axes_length(:,i),mymovie.(type).orientations(1,i));
          angl = ell_pos(cyt1,1);
          dist = abs(ell_pos(:,1)-(2*pi - angl));
          second_indx = find(dist == min(dist) & dist < angl_thresh);

          if (~isempty(second_indx))
            cyt2 = mymovie.(type).ruffles(i).cluster(second_indx,1);
          end
        end
      end

      if (cyt1 ~= 0)
        cyt1 = mymovie.(type).ruffles(i).cluster(cyt1,1);
      end
      %if (cyt1 == 0 & cyt2 == 0)
      %  val1 = 0;
      %  val2 = 0;
      %end
    end
    
    %values(:,i) = [val1;val2];

    %if ((val1 < thresh && val2 < min_thresh) || (val2 < thresh && val1 < min_thresh))
    if (cyt1 == 0 & cyt2 == 0)
      cytokinesis = i;
      disp(['Cytokinesis at ' num2str(i)]);
      mymovie.(type).cytokinesis = cytokinesis;
      times = frame_timing(mymovie);

      if (~isempty(times))
        mymovie = frame_timing(mymovie, -times(cytokinesis));
      end

      break;
    end
  end

  %figure;
  %plot(a)
  %hold on;
  %plot(values(1,:),'r')
  %plot(values(2,:),'g')

  %keyboard

  return;
end
