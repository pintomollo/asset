function [dist, weight] = spot_similarity(pts1, pts2)
      
  dist = sqrt(bsxfun(@minus,pts1(:,1),pts2(:,1).').^2 + bsxfun(@minus,pts1(:,2),pts2(:,2).').^2);

  if (size(pts1,2) > 2)
    signal1 = pts1(:,3).^2 .* pts1(:,4);
    signal2 = pts2(:,3).^2 .* pts2(:,4);

    weight = max(bsxfun(@rdivide,signal1,signal2.'), bsxfun(@rdivide,signal1,signal2.')).^2;
  else
    weight = ones(size(dist));
  end

  return;
end
