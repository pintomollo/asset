function [params] = fit_front(func, p0, x, y, lbound, ubound, opts)

  if (exist('fit_front_mex') == 3)
    y_tmp = y;
    MaxIter = 400;
    stop_tol = 1e-6;
    [niter, params] = fit_front_mex(p0, y_tmp, MaxIter, stop_tol, lbound, ubound, x);
  else
    params = lsqcurvefit(func, p0, x, y, lbound, ubound, opts);
  end

  return;
end
