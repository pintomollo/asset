function [clusters] = cluster_distance(positions, dist_thresh, max_clusts)

  dists = pdist(positions);
  links = linkage(dists, 'average');
  clusters = cluster(links, 'cutoff', dist_thresh, 'criterion', 'distance');

  vals = sort(clusters);
  [clust, indexes] = unique(vals);
  nelems = diff([indexes; length(vals)]);
  [junk, mapping] = sort(nelems, 'descend');
  [junk, mapping] = sort(mapping);

  clusters = mapping(clusters);

  return;
end
