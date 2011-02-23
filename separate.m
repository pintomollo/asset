function img = separate(img)
  %figure;imshow(img);
  global USE_MM;

  if (~any(any(img)))
    return;
  end

  if (USE_MM)
    dist = double(mmdist(img));
  else
    dist = bwdist(~img);
  end

  dist = -dist;
  dist(~img) = -Inf;

  pixs = sort(dist(dist ~= -Inf).');
  thresh = pixs(floor(length(pixs)/20)+1);
  dist(dist<thresh) = -Inf;

  %figure;imshow(dist,[],'InitialMagnification','fit');

  if (USE_MM)
    dist = dist - thresh + 1;
    dist(isinf(dist)) = 0;
    labels = mmwatershed(uint16(dist), mmsebox, 'REGIONS') - 1;
  else
    labels = watershed(dist) - 1;
  end 

  %figure;imshow(labels,[],'InitialMagnification','fit');
  prohib = unique(labels(~img));

  narea = max(max(labels));
  sizes = zeros(1, narea);

  for i=1:narea
    if (~any(i==prohib))
      sizes(i) = sum(sum(labels==i));
    end
  end
  
  winner = find(sizes==max(sizes),1);

  if (isempty(winner))
    return;
  end

  img = (labels==winner);

  %figure;imshow(img);

end


