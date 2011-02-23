function best = fit_logistic(x, y, x0)

  
    opts = optimset('Display','off', 'Algorithm', 'levenberg-marquardt', 'MaxFunEvals', 5000, 'MaxIter', 5000);
    [best] = lsqnonlin(@fit_error, x0, [], [], opts);

  return;

  function err = fit_error(p)
    
    t = logistic(x, p(1), p(2), p(3), p(4), p(5));
    err = sum(abs(t - y));

    return;
  end
end
