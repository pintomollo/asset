function mymovie = align_embryo(mymovie, opts)

  if (isfield(mymovie, 'markers') & ~isempty(mymovie.markers) & isfield(mymovie.markers, 'orientations') & ~isempty(mymovie.markers.orientations))
    type = 'markers';
    if (isfield(mymovie, 'eggshell') & ~isempty(mymovie.eggshell))
      [nframes imgsize] = size_data(mymovie.eggshell);
    else
      [nframes imgsize] = size_data(mymovie.cortex);
    end
  elseif (isfield(mymovie, 'dic') & ~isempty(mymovie.dic) & isfield(mymovie.dic, 'orientations') & ~isempty(mymovie.dic.orientations))
    type = 'dic';
    [nframes, imgsize] = size_data(mymovie.(type));
  elseif (isfield(mymovie, 'data') & ~isempty(mymovie.data) & isfield(mymovie.data, 'orientations') & ~isempty(mymovie.data.orientations))
    type = 'data';
    [nframes, imgsize] = size_data(mymovie.(type));
  else
    error 'None of the expected field are present in ''mymovie''';
  end
  opts.segmentation_type = type;

  mymovie.(type) = align_orientations(mymovie.(type));
  mymovie = find_ruffles(mymovie, opts);

  if (isfield(mymovie.(type), 'inverted'))
    actu_invert = mymovie.(type).inverted;
  else
    actu_invert = false;
  end

  nbins = 16;

  pts = zeros(0, 3);

  for i=1:nframes
    carth = mymovie.(type).ruffles(i).carth;

    if (isempty(carth))
      continue;
    end

    %if (i==50)
    %  keyboard
    %end
    ell = carth2elliptic(carth, mymovie.(type).centers(:,i),mymovie.(type).axes_length(:,i),mymovie.(type).orientations(1,i));
    ell = [ell i*ones(size(ell, 1), 1)];
    pts = [pts; ell];
  end

  bins = [0:2*pi/nbins:2*pi + 1e-5];
  [counts, indx] = histc(pts(:,1), bins);

  myhist = counts(1:nbins/2) + counts(end-1:-1:end-nbins/2);

  inverted = false(3,1);
  inverted(1) = (sum(myhist(1:nbins/8)) > sum(myhist(end-(nbins/8)+1:end)));

  [junk, cytok_pos] = max(myhist(nbins/8+1:end-(nbins/8)));
  inverted(2) = (cytok_pos > nbins/8);

  if (inverted(1) ~= inverted(2))
    angl_thresh = 0.2;

    for i=nframes:-1:1
    
      ell_pos = pts(pts(:,3)==i, 1:2);

      if (isempty(ell_pos))
        continue;
      end

      dists = abs(ell_pos(:,1) - pi/2 - pi*(ell_pos(:,1) > pi));
      max_indx = find(dists == min(dists) & dists < angl_thresh, 1);
      if (isempty(max_indx))
          continue;
      end
      
      angl = ell_pos(max_indx, 1); 

      inverted(3) = ((angl < pi & angl > pi/2) | (angl > pi & angl < 3*pi/2));
      break;
    end
  end

  orient = mymean(mymovie.(type).orientations);
  if (sum(inverted) > 1)
    if (opts.verbosity > 0)
      disp(['A-P were inverted !']);
    end
    inverted = ~actu_invert;

    orient = orient + pi;
    orient = orient - 2*pi*(orient > 2*pi);

    mymovie.(type) = align_orientations(mymovie.(type), orient);
  else
    inverted = actu_invert;
  end

  fields = fieldnames(mymovie);
  for f = 1:length(fields)
    if (~isempty(mymovie.(fields{f})) && isfield(mymovie.(fields{f}), 'orientations') && ~isempty(mymovie.(fields{f}).orientations))

      mymovie.(fields{f}) = align_orientations(mymovie.(fields{f}), orient);
      mymovie.(fields{f}).inverted = inverted;

      %if ((~isfield(mymovie.(fields{f}), 'inverted') & inverted) | (isfield(mymovie.(fields{f}), 'inverted') & mymovie.(fields{f}).inverted ~= inverted))
      %  tmp_orient = mymovie.(fields{f}).orientations;
      %  tmp_orient = tmp_orient + pi;
      %  tmp_orient = tmp_orient - 2*pi*(tmp_orient > 2*pi);
      %  mymovie.(fields{f}).orientations = tmp_orient;
      %  mymovie.(fields{f}).inverted = inverted;
      %end
    end
  end

  if (isfield(mymovie, 'metadata') & ~isempty(mymovie.metadata) & isfield(mymovie.metadata, 'orientation_3d') & ~isempty(mymovie.metadata.orientation_3d))
    mymovie.metadata.orientation_3d = align_orientations(mymovie.metadata.orientation_3d, orient);
  end


  %depth = zeros(size(counts));
  %for i=1:length(bins)
  %  depth(i) = sum(pts(indx == i, 2));
  %end

  %figure;scatter(bins, counts);
  %hold on;
  %scatter(bins, depth, 'r');
  %keyboard

  return;
end
