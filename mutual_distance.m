function dist = mutual_distance(pts1, pts2)
      
  dist = sqrt(bsxfun(@minus,pts1(:,1),pts2(:,1).').^2 + bsxfun(@minus,pts1(:,2),pts2(:,2).').^2);

  return;
end
