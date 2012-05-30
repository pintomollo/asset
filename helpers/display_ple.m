function display_ple(chi2ple, psple)

  alpha = 0.68;

  nparams = length(chi2ple);
  chi1 = chi2inv(alpha, 1);
  chin = chi2inv(alpha, nparams);

  figure;
  for i=1:nparams
    pos = (psple{i}(:, i));

    subplot(1, nparams, i);
    semilogx(pos, chi2ple{i});
    hold on;
    semilogx(pos,chi1*ones(size(pos)));
    semilogx(pos,chin*ones(size(pos)));
  end

  return;
end
