function img = separate(img)

  if (~any(any(img)))
    return;
  end

  dist = bwdist(~img);

  dist = -dist;
  dist(~img) = -Inf;

  pixs = sort(dist(dist ~= -Inf).');
  thresh = pixs(floor(length(pixs)/20)+1);
  dist(dist<thresh) = -Inf;

  %figure;imshow(dist,[],'InitialMagnification','fit');

  labels = watershed(dist) - 1;

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


