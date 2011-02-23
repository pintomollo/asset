function [dist, weight, alt_weight] = merging_distance(pts1, pts2, spots)

  dist = sqrt(bsxfun(@minus,pts1(:,1),pts2(:,1).').^2 + bsxfun(@minus,pts1(:,2),pts2(:,2).').^2);
  weight = ones(size(dist));

  alt_weight = weight.';

  return;
end
