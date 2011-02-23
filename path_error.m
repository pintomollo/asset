function errors = path_error(values, paths, nbins)

  [npts, nsplines, ngroups] = size(paths);
  dt = diff(values);

  theta = (2*pi) / nbins;
  bins = [0 0.5:nbins nbins] * theta;

  for g = 1:ngroups
    for s = 1:nsplines
      dr = paths(:,s,g);
      ds = abs(diff(sign(dr)));
      dr = abs(dr);

      areas = (dr(1:end-1) + dr(2:end)) .* dt ./ 2;

      if (any(ds ~= 0))
        targets = (ds ~= 0);
        denum = dt(targets) ./ (2* (dr([targets; false]) + dr([false; targets])));
        indxs = ~isinf(denum);

        denum = denum(indxs);
        targets(~indxs) = false;
        
        areas(targets) = (dr([targets; false]).^2 + dr([false; targets]).^2) .* denum;
      end

      [junk, map] = histc(values, bins);

      for k=1:nbins
        errors(g,s,k) = sum(areas(map == k));
      end
      errors(g,s,1) = errors(g,s,1) + sum(areas(map(1:end-1) == (nbins + 1)));

    end
  end

  errors = errors / (2*pi);

  return;
end
