function display_mcmc_results(fname)

  [header, data, varying] = parse_mcmc_results(fname);
  nparams = sum(varying);

  if (size(header, 2) == 4)
    valids = ~(header(:, 4))
    header = header(valids, :);
    data = data(valids, :);
  end

  return;
end
