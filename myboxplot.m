function myboxplot(varargin)

  pts = [];
  groups = [];

  for i=1:nargin
    vals = varargin{i}(:);
    pts = [pts; vals];
    groups = [groups; i*ones(size(vals))];
  end
  figure;boxplot(pts, groups);

  [h,p] = myttest(pts, groups)

  return;
end
