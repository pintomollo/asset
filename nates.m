function du = nates(u, h, params)

  int_trapz = (h/2) * (sum(u, 2) + sum(u(:, 2:end-1), 2));
  cyto = bsxfun(@minus, params([5 10]).', params(11)*int_trapz/params(12));

  du = bsxfun(@minus, bsxfun(@times, params([1 6]).', cyto), bsxfun(@times, params([2 7]).', u)) - ...
       bsxfun(@times, params([3 8]).', bsxfun(@power, u([2 1], :), params([4 9]).') .* u);

  return;
end
