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

  a = mymean(pts, 1, groups)

  a(2) / a(1)

  mean(a(3:4) ./ a(1:2))
  mean(a(5:6) ./ a(1:2))

  return;
end
