% starting from optimized (the best you could get) parameter values <p>,
% with current value of the objective function <chi2>,
% calculate the profile likelihood for parameter with index <jk>,
% taking <samplesize> steps up and down.
% The function returns <chi2ple>, the values of the profile likelihood, and <psple>,
% the corresponding parameter values. 

function [chi2ple, psple] = ple(p, jk, samplesize, func, varargin)

  % Initialize the likelihood for the current position
  chi2 = func(p, varargin{:});

  % remember current values
  pbackup = p;
  chi2backup = chi2;

  % degrees of freedom desired for confidence level
  %df = 1; % corresponds to point-wise confidence intervals
  df = length(p); % corresponds to simultaneous confidence intervals

  % desired level of confidence
  alpha = 0.95;

  % <alpha> quantile of chi2-distribution with <df> degrees of freedom 
  chi2_threshold = chi2inv(alpha, df);

  do_all = true;
  step_thresh = -1;
  if (~isfinite(samplesize))
    do_all = false;
    step_thresh = (chi2_threshold - chi2) * 1e-10;

    % Some fairly large number !
    samplesize = 1e6;
  end

  % initialize vector to store the values of the objective function
  % along the profile likelihood
  chi2ple = NaN(2 * samplesize + 1, 1);

  % initialize matrix to store the parameter values along the profile likelihood
  psple = NaN(2 * samplesize + 1, length(p));

  % set current parameter values and value of the objective function
  % as central values
  chi2ple(samplesize + 1) = chi2;
  psple(samplesize + 1, :) = p;

  % calculate profile likelihood for increasing value of p(jk)
  dpjk = p(jk)*0.1; % inital_step() increases p(jk) by 10% initially
  for i = 1:samplesize
    % increase p(jk)
    [p, dpjk] = init_step_direct(chi2, p, jk, dpjk, chi2_threshold, func, varargin{:});

    % fit parameters EXCEPT for parameter with index jk
    [p, chi2] = fit_func(p, jk, func, varargin{:});

    % store values
    chi2ple(samplesize + 1 + i) = chi2;
    psple(samplesize + 1 + i, :) = p;

    fprintf(1, '+');

    if (~do_all & (abs(chi2ple(samplesize + 1 + i) - chi2ple(samplesize + i)) < step_thresh | chi2 > chi2_threshold))
      break;
    end
  end

  % reset original values
  p = pbackup;
  chi2 = chi2backup;

  % calculate profile likelihood for decreasing value of p(jk)
  dpjk = -p(jk)*0.1; % inital_step() decreases p(jk) by 10% initially
  for i = 1:samplesize
    % decrease p(jk)
    p = init_step_direct(chi2, p, jk, dpjk, chi2_threshold, func, varargin{:});

    % fit parameters EXCEPT for parameter with index jk
    [p, chi2] = fit_func(p, jk, func, varargin{:});

    % store values
    chi2ple(samplesize + 1 - i) = chi2;
    psple(samplesize + 1 - i, :) = p;
    fprintf(1, '-');

    if (~do_all & (abs(chi2ple(samplesize + 1 - i) - chi2ple(samplesize + 2 - i)) < step_thresh | chi2 > chi2_threshold))
      break;
    end
  end
  fprintf(1, '\n');

  if (~do_all)
    valids = ~isnan(chi2ple);
    chi2ple = chi2ple(valids);
    psple = psple(valids, :);
  end

  return;
end

% increase/decrease value of p(jk)
function [p, dpjk] = init_step_direct(chi2, p, jk, dpjk, chi2_threshold, func, varargin)
  % percentage increase of objective function allowed for the step
  chi2_relative_increase = 0.1;

  % increase of objective function due to step
  ptmp = p;
  ptmp(jk) = p(jk) + dpjk;
  chi2_diff = func(ptmp, varargin{:}) - chi2; % evaluate the objective function

  step_min = chi2_diff / 1e4;
  prev_step = 0;
  count_max = 100;
  count = 0;

  if(chi2_diff > chi2_threshold * chi2_relative_increase) % decrease dpjk
    while(chi2_diff > chi2_threshold * chi2_relative_increase & count < count_max & abs(chi2_diff - prev_step) > step_min)
      dpjk = dpjk / 2; % halve step increment
      ptmp(jk) = p(jk) + dpjk;
      prev_step = chi2_diff;
      chi2_diff = func(ptmp, varargin{:}) - chi2; % evaluate the objective function

      count = count + 1;
    end
  else % increase dpjk
    while(chi2_diff < chi2_threshold * chi2_relative_increase & count < count_max & abs(chi2_diff - prev_step) > step_min)
      dpjk = dpjk * 2; % double step increment
      ptmp(jk) = p(jk) + dpjk;
      prev_step = chi2_diff;
      chi2_diff = func(ptmp, varargin{:}) - chi2; % evaluate the objective function

      count = count + 1;
    end
    dpjk = dpjk / 2; % halve step increment
  end

  % finally set new parameter values
  p(jk) = p(jk) + dpjk;
end

% fitting interface to your own code
% for fitting parameters <p> EXCEPT for parameter with index <jk>
function [p, chi2] = fit_func(p, jk, func, varargin)

  opts = optimset('Display','off', 'Algorithm', 'levenberg-marquardt', 'MaxFunEvals', 5000, 'MaxIter', 5000);

  index = [1:length(p)];
  fitted = (index ~= jk);
  [best, junk, chi2] = lsqnonlin(@generic_func, p(fitted), [], [], opts);
  p(fitted) = best;

  return;

  function val = generic_func(p_tmp)

    new_p = p;
    new_p(fitted) = p_tmp;
    val = func(new_p, varargin{:});

    return;
  end
end
