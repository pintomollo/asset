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
  opts.recompute = false;
  mymovie = find_ruffles(mymovie, opts);

  if (isfield(mymovie.(type), 'inverted'))
    actu_invert = mymovie.(type).inverted;
  else
    actu_invert = false;
  end

  nbins = 32;

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

    if (isempty(ell))
      continue;
    end
    ell = [ell(:,1) mymovie.(type).ruffles(i).properties(:,1) i*ones(size(ell, 1), 1)];
    pts = [pts; ell];
  end

  if (isempty(pts))
    warning('Orientation of the embryo cannot be determined due to the lack of ruffles.')

    return;
  end

  %figure;scatter3(pts(:,1), pts(:,2), pts(:,3));

  bins = [0:2*pi/nbins:2*pi + 1e-5];
  [counts, indx] = histc(pts(:,1), bins);

  valids = (indx > 0);
  depths = zeros(nbins, 1);
  depths(counts > 0) = mymean(pts(valids,2), 1, indx(valids));

  myhist = counts(1:nbins/2) + counts(end-1:-1:end-nbins/2);
  mydepth = depths(1:nbins/2) + depths(end-1:-1:end-nbins/2);
  mydepth2 = cumsum(mydepth(end:-1:1)); 

  inverted = false(3,1);

  %figure;subplot(121);plot(myhist);subplot(122);plot(mydepth)

  %% Couting invaginations on both sides
  inverted(1) = (sum(myhist(1:nbins/8)) > sum(myhist(end-(nbins/8)+1:end)));
  %inverted(1) = (sum(mydepth(1:nbins/8)) > sum(mydepth(end-(nbins/8)+1:end)));

  %% Checking where the deepest ones are located
  % And the smalle ones, compare basically smooth vs rough
  mydepth = cumsum(mydepth);
  mydepth = (mydepth / mydepth(end));
  inverted(2) = (mydepth(nbins/4) < 0.5);

  %[junk, cytok_pos] = max(mydepth(nbins/8+1:end-(nbins/8)))
%  [junk, cytok_pos] = max(mydepth);
%  inverted(2) = (cytok_pos >= nbins/4);

%  if (inverted(1) ~= inverted(2))

    times = get_manual_timing(mymovie, opts);
    frames = [];
    if (~isnan(times(3)))
      frames = [frames times(3):nframes];
    elseif (~isnan(times(1)))
      frames = [frames times(1)-5:times(1)+5];
    else
      frames = [1:nframes];
    end
    
    frames = frames(frames > 0 & frames <= nframes);

    good_ruffles = (ismember(pts(:,3), frames) & pts(:,2)>0.1);
    if(~any(good_ruffles))
      good_ruffles = (ismember(pts(:,3), frames));
    end

    if (any(good_ruffles))
      targets = pts(good_ruffles, :);
      targets(targets(:,1)>pi,1) = 2*pi - targets(targets(:,1)>pi,1);

      weights = targets(:,2) .* (targets(:,3) / frames(end));
      weights = weights / sum(weights);
      %avg_angl = sum(targets(:,1) .* weights);

      % Weighted median
      [angls, indxs] = sort(targets(:,1));
      weights = cumsum(weights(indxs));
      indx = find(weights>0.5, 1, 'first');
      med_angl = angls(indx);

      inverted(3) = (med_angl > pi/2);
    end

  %keyboard

    %{
    angl_thresh = 0.2;
    times = get_manual_timing(mymovie, opts);
    tests = NaN(nframes, 1);
    angls = NaN(nframes, 1);
    max_depth = 0;

    for i=nframes:-1:1
    
      ell_pos = pts(pts(:,3)==i, 1:2);
      npos = size(ell_pos, 1);

      if (npos == 0)
        continue;
      end

      dists = abs(ell_pos(:,1) - pi/2 - pi*(ell_pos(:,1) > pi));
      %[dists, dindxs] = sort(dists);
      [depths, pindxs] = sort(ell_pos(:,2), 1, 'descend');
      %[junk, pindxs] = sort(pindxs);
      %inv_pindxs = pindxs(dindxs);

      %if (inv_pindxs(1) > 2)
      %  continue;
      %end

      %curr_depth = depths(inv_pindxs(1));
      curr_depth = depths(1);
      if (curr_depth > max_depth)
        max_depth = curr_depth;
      end

      angl = ell_pos(pindxs(1), 1);
      if (numel(dists) > 1)
        %if (inv_pindxs(1) == 1)
          if (depths(2) > depths(1)/2)
            angl = ((2*pi - ell_pos(pindxs(2), 1)) + angl) / 2;
          end
        %else
        %  if (depths(1) > depths(2)/2)
        %    angl = ((2*pi - ell_pos(pindxs(1), 1)) + angl) / 2;
        %  end
        %end
      end
      tests(i) = ((angl < pi & angl > pi/2) | (angl > pi & angl < 3*pi/2));
      angls(i) = angl;

%      max_indx = find(dists == min(dists) & dists < angl_thresh, 2)
%      if (isempty(max_indx))
%          continue;
%      end
      
%      angl = ell_pos(max_indx, 1); 
%      [angl ell_pos(max_indx, 2)]

%      inverted(3) = ((angl < pi & angl > pi/2) | (angl > pi & angl < 3*pi/2));
%      break;
      if (curr_depth < max_depth/2 | i < times(end))
      %if (i < times(end))
        break;
      end
    end
    inverted(3) = (mymean(tests)+(1e-6*(tests(end)-0.5)) > 0.5);
    %}

    %[junk, smooth_pos] = min(mydepth);
    %inverted(3) = (smooth_pos >= nbins/4);
%  end

%  if (any(inverted))
%    beep;
%    keyboard
%  end

  orient = mymean(mymovie.(type).orientations);
  if (sum(inverted) > 1)
    if (opts.verbosity > 0)
      disp(['A-P were inverted !']);
    end

    accept_detection = true;
    if (actu_invert)
      if (opts.verbosity > 2)
        answer = questdlg('The orientation of this embryo was previously set manually, revert to the detected one instead ?', 'ASSET');
      elseif (opts.verbosity > 0)
        answer = input('The orientation of this embryo was previously set manually,\nrevert to the detected one instead ? Y/[N]:', 's');
      else
        answer = 'no';
      end

      accept_detection = strncmpi(answer, 'Y', 1);
    end

    if (accept_detection)
      inverted = false;

      orient = orient + pi;
      orient = orient - 2*pi*(orient > 2*pi);
    else
      inverted = actu_invert;
    end

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
