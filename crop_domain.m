function [domain, boundary] = crop_domain(domain, center)

  boundary = [];
  goods = ~isnan(domain);

  for i=1:size(domain, 1);
    if (~any(goods(i,:)))
      goods(i,:) = true;
    else
      alls = find(goods(i,:));
      first = alls(1);
      last = alls(end);

      goods(i,first:last) = true;
    end
  end
  valids = all(goods, 1);

  alls = find(valids);
  first = alls(1);
  last = alls(end);

  if (nargin == 2)
    first = center - first;
    last = last - center;

    boundary = min(first, last);
    domain = domain(:, [-boundary:boundary]+center);
  else
    domain = domain(:, first:last);
  end

  return;
end
