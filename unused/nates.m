function du = nates(u, h, params, fid)

  int_trapz = (h/2) * (sum(u, 1) + sum(u(2:end-1, :), 1));
  cyto = bsxfun(@minus, params(5, :), params(end-1)*int_trapz/params(end));

  du = bsxfun(@minus, bsxfun(@times, params(1, :), cyto), bsxfun(@times, params(2,:), u)) - ...
       bsxfun(@times, params(3,: ), bsxfun(@power, u(:, [2 1]), params(4, :)) .* u);

%  if (nargin == 4)
%    fprintf(fid, '%f=%f=%f=%f=%f=%f=', sum(u(1,2:end-1), 2), sum(u(1,:), 2), int_trapz(1), sum(u(2,2:end-1), 2), sum(u(2,:), 2), int_trapz(2));
%    fprintf(fid, ':(%f,%f,%f):', params(5), int_trapz(1), params(11)/params(12));
%  fprintf(fid, '%f %f\n', du);
%      fprintf(fid, '-----------\n');
%  end

  return;
end
