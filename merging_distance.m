function [dist, weight, alt_weight] = merging_distance(pts1, pts2, spots, links)

  dist = sqrt(bsxfun(@minus,pts1(:,1),pts2(:,1).').^2 + bsxfun(@minus,pts1(:,2),pts2(:,2).').^2);

  [nends, ncols] = size(pts1);
  ninterm = size(pts2, 1);

  int_indx = 3;
  %int_indx = min(4, ncols-2);

  frames = unique(pts2(:, end)).';
  prev_intense = NaN(ninterm, 1);
  for i = frames
    tmp_links = links{i};
    tmp_vals = NaN(max(tmp_links(:,1)), 1);
    tmp_vals(tmp_links(:,1)) = spots{i-1}(tmp_links(:,2), int_indx);

    currents = (pts2(:,end) == i);
    prev_intense(currents) = tmp_vals(pts2(currents, end-1));
  end

  weight = bsxfun(@rdivide, bsxfun(@plus,pts1(:,int_indx),prev_intense(:,1).'), pts2(:, int_indx).');
  weight(weight < 1) = weight(weight < 1).^(-2);

  alt_weight = pts2(:, int_indx) ./ prev_intense;
  alt_weight(alt_weight < 1) = alt_weight(alt_weight < 1).^(-2);
  alt_weight(isnan(alt_weight)) = max(alt_weight, [], 1);
  alt_weight = diag(alt_weight);
  alt_weight(alt_weight == 0) = Inf;

  return;
end
