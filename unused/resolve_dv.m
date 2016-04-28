function mymovie = resolve_dv(mymovie)

  nframes = length(mymovie.data.centrosomes);

  thresh = 5;
  pos = NaN(nframes,2);

  mymovie = carth2normalized(mymovie);

  for i=1:nframes
    if (any(~isnan(mymovie.data.centrosomes(i).carth(:))))
      pos(i,:) = mymovie.data.centrosomes(i).warped(2,:);
    end
  end

  pos = pos(~isnan(pos(:,1)),:);

  if (isempty(pos))
    return;
  end

  if (size(pos,1) > 6)
    init_pos = mean(pos(1:3,:),1);
    end_pos = mean(pos(end-2:end,:),1);

    other = pos(4:end-3,:);
    dinit = sqrt((other(:,1)-init_pos(1)).^2 + (other(:,2)-init_pos(2)).^2);
    dend = sqrt((other(:,1)-end_pos(1)).^2 + (other(:,2)-end_pos(2)).^2);
    
    pos = pos(dinit > thresh & dend > 5, :);
  end

  above = sum(pos(:,2) < 0) / size(pos, 1);

  if (above < 0.5)
    disp('D-V were inverted !');
    mymovie.data.dv_inverted = ~mymovie.data.dv_inverted;
    mymovie = carth2normalized(mymovie);
  %else
  %  mymovie.data.dv_inverted = false;
  end

  return;
end
