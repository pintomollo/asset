function sample_size = estimate_sample_size(ci, precision)
% [Greenwood & Sandomire, Sample Size Required for Estimating the Standard Deviation as a Per Cent of its True Value, JASA, 1950]

  n = 10.^[0:4];

  p_plus = (1+precision)^2;
  p_minus = (1-precision)^2;

  ps = curr_prob(n);

  cross_pos = (abs(diff(sign(ps))) > 0); 
  if (length(cross_pos) > 0 & any(cross_pos))
    indx = find(cross_pos, 1);
    ns = n([indx indx+1]);
  else
    if (all(ps < 0 | ~isfinite(ps)))
      sample_size = Inf;
    else
      sample_size = 1;
    end

    return;
  end

  N = fzero(@curr_prob, ns);

  sample_size = ceil(N) + 1;

  return;

  function prob_precision = curr_prob(dof)

    chi_high = dof*p_plus;
    chi_low = dof*p_minus;

    p1 = chi2cdf(chi_high, dof);
    p2 = chi2cdf(chi_low, dof);

    prob_precision = ci - (p1 - p2);

    return;
  end
end
