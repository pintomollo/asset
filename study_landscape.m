function study_landscape(file)

  %if (ischar(file))
  %  file = load(file);
  %end

  fitting = get_struct('fitting');
  fitting = load_parameters(fitting, 'param_distribution.txt');

  init_pos = fitting.init_pos;
  variation = [0.66 1.5];
  nparams = length(init_pos);

  for i = 1:nparams
    for v = variation
      tmp_pos = init_pos;
      tmp_pos(i) = tmp_pos(i)*v;
      find_kymograph(file, 'fit_full', true,  'config_fitting', 'param_distribution', 'integrate_sigma', true, 'init_pos', tmp_pos, 'max_iter', 2000);
    end
  end

  return;
end
