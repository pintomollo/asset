function [fit_params, fit_energy, fit_temperatures, fit_viscosity] = model_description(parameter_set)

  fit_energy = [];
  fit_temperatures = false;
  fit_viscosity = false;

  switch parameter_set
    case 0
      fit_params = [];
    case 1
      fit_params = [4 12];
    case 2
      fit_params = [4 5 12 13];
    case 12
      fit_params = [4 5 12 13];
      fit_temperatures = true;
      fitting.fit_model = true;
    case 13
      fit_params = [4 5 12 13];
      fit_temperatures = true;
      fit_energy = 0.65;
    case 14
      fit_params = [4 5 12 13];
      fit_temperatures = true;
      fit_energy = [0.65 0.65];
    case 15
      fit_params = [4 5 12 13];
      fit_temperatures = true;
      fit_energy = [0.65 0.65 0.84];
    case 16
      fit_params = [4 5 12 13];
      fit_temperatures = true;
      fit_energy = ones(1,4) * 0.65;
    case 17
      fit_params = [4 5 12 13];
      fit_temperatures = true;
      fit_energy = ones(1,7) * 0.65;
    case 20
      fit_params = [4 5 12 13];
      fit_viscosity = true;
      fit_temperatures = false;
      fitting.fit_model = true;
    case 22
      fit_params = [4 5 12 13];
      fit_viscosity = true;
      fit_temperatures = true;
      fitting.fit_model = true;
    case 23
      fit_params = [4 5 12 13];
      fit_viscosity = true;
      fit_temperatures = true;
      fit_energy = 0.65;
    case 24
      fit_params = [4 5 12 13];
      fit_viscosity = true;
      fit_temperatures = true;
      fit_energy = [0.65 0.65];
    case 25
      fit_params = [4 5 12 13];
      fit_viscosity = true;
      fit_temperatures = true;
      fit_energy = ones(1,4) * 0.65;
    case 26
      fit_params = [4 5 12 13];
      fit_viscosity = true;
      fit_temperatures = true;
      fit_energy = ones(1,7) * 0.65;
    case 30
      fit_params = [1:14];
      fit_temperatures = true;
      fit_energy = [0.65 0.65 0.84];
      %fit_viscosity = true;
      %fit_energy = [0.65 0.65];
    case 3
      fit_params = [4 5 6 12 13];
    case 4
      fit_params = [2 4 5 6 10 12 13];
    case 5
      fit_params = [4 5 6 12 13 14];
    case 6
      fit_params = [1:14];
    case 7
      fit_params = [1:6 9:13];
    otherwise
      fit_params = [1];
  end

  return;
end
