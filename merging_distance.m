function [dist, weight, alt_weight] = merging_distance(pts1, pts2, spots, links)

  dist = sqrt(bsxfun(@minus,pts1(:,1),pts2(:,1).').^2 + bsxfun(@minus,pts1(:,2),pts2(:,2).').^2);

  nends = size(pts1, 1);
  ninterm = size(pts2, 1);

  frames = unique(pts2(:, end)).';
  prev_intense = NaN(ninterm, 1);
  for i = frames
    tmp_links = links{i};
    tmp_vals = NaN(max(tmp_links(:,1)), 1);
    tmp_vals(tmp_links(:,1)) = spots{i-1}(tmp_links(:,2), 4);

    currents = (pts2(:,end) == i);
    prev_intense(currents) = tmp_vals(pts2(currents, end-1));
  end

  weight = bsxfun(@rdivide, bsxfun(@plus,pts1(:,4),prev_intense(:,1).'), pts2(:, 4).');
  weight(weight < 1) = weight(weight < 1).^(-2);

  alt_weight = pts2(:, 4) ./ prev_intense;
  alt_weight(alt_weight < 1) = alt_weight(alt_weight < 1).^(-2);
  alt_weight = diag(alt_weight);
  alt_weight(alt_weight == 0) = Inf;

  return;
end
