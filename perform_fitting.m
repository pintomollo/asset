function perform_fitting(selection)

  params = {'config_fitting'; 'fit_kymo'; 'config_modeling'; 'custom_flow'};

  switch selection
    case 0
      % All ani-2 and c27d9.1
      files = [dir('1056-3-*_.mat'); ...
               dir('1056-14-*_.mat'); ...
               dir('1056-24-*_.mat'); ...
               dir('1056-ani2-*_.mat'); ...
               dir('1056-c27d91-*_.mat')];
      repeats = 1;
      init_noise = 0;
      starts = 'B';
      param_set = 2;
    case 1
      % All 24°C
      files = dir('1056-24-*_.mat');
      repeats = 1;
      init_noise = 0.1;
      starts = 'GABB';
      param_set = 2;
    case 2
      % All 14°C
      files = dir('1056-14-*_.mat');
      repeats = 1;
      init_noise = 0.1;
      starts = 'GABB';
      param_set = 2;
    case 3
      % All 3°C, ani-2 and c27d9.1
      files = [dir('1056-3-*_.mat'); ...
               dir('1056-ani2-*_.mat'); ...
               dir('1056-c27d91-*_.mat')];
      repeats = 1;
      init_noise = 0.1;
      starts = 'GABB';
      param_set = 2;
    case 4
      files = dir('1056-*-all.mat');
      repeats = 5;
      init_noise = 0.1;
      starts = 'GABB';
      param_set = 2;
    case 5
      files = {'1056-temps-all.mat'};
      repeats = 2;
      init_noise = 0.1;
      starts = 'GABB';
      param_set = [12 13 23];
    case 6
      files = {'1056-temps-all.mat'};
      repeats = 2;
      init_noise = 0.1;
      starts = 'GABB';
      param_set = [24 25];
    case 7
      files = {'simulation'};
      repeats = 10;
      init_noise = [0.1 0.25 0.5 1];
      starts = 'GA';
      param_set = 2;
    case 8
      perform_fitting(8.1);
      perform_fitting(8.2);
      return;

    case 8.1
      files = dir('1056-*-all.mat');
      repeats = 1;
      init_noise = 0;
      starts = 'O';
      param_set = 6;
      params{2} = 'hessian';
    case 8.2
      files = dir('1056-*-all.mat');
      repeats = 1;
      init_noise = 0;
      starts = 'O';
      param_set = 6;
      params{2} = 'sensitivity';
    case 9
      files = {'1056-temps-all.mat'};
      repeats = 1;
      init_noise = 0;
      starts = 'T';
      param_set = [12 13 23 24 25];
      params{2} = 'fit_temperature';
    case 10
      files = {'1056-all-all.mat'};
      repeats = 2;
      init_noise = 0;
      starts = 'T';
      param_set = 2;
    case 11
      files = {'1056-all-all.mat'};
      repeats = 2;
      init_noise = 0;
      starts = 'T';
      param_set = 13;
    case 12
      files = {'1056-all-all.mat'};
      repeats = 2;
      init_noise = 0;
      starts = 'T';
      param_set = 24;
    case 13
      files = {'1056-temps-all.mat'};
      repeats = 2;
      init_noise = 0;
      starts = 'T';
      param_set = 24;
      params{2} = 'fit_flows';
    case 14
      files = {'1056-temps-all.mat'};
      repeats = 2;
      init_noise = 0.1;
      starts = 'T';
      param_set = 24;
      params{2} = 'fit_temp_indep';

    otherwise
      warning('Choose a fitting group between 1 and 7');

      return;
  end

  ntotal = length(files)*repeats*length(init_noise)*length(starts)*length(param_set);
  counts = 1;

  for s = 1:length(starts)
    switch starts(s)
      case 'A'
      % Best average 20°C detected previously
        s_params = {'init_pos'; [0.00107 1.9688 0.0069 2.4249]};
      case 'B'
      % Best value found
        s_params = {'start_with_best'; true};
      case {'O', 'T'}
      % Current best found during optimization in BestFits
        s_params = {'init_pos'; []};
      otherwise
      % Default setting
        s_params = {};
    end

    for i = 1:length(init_noise)
      i_params = {'init_noise'; init_noise(i)};

      for p = 1:length(param_set)
        p_params = {'parameter_set'; param_set(p)};

        for r = 1:repeats

          for f = 1:length(files)
            if (iscell(files))
              f_params = files(f);
            else
              f_params = {files(f).name};
            end

            if (strncmp(f_params{1}, 'simulation', 10))
              f_params = [f_params; 'type'; f_params];
            end

            if (starts(s) == 'O')
              vals = group_ml_results('BestFits/adr-kymo-*_evol.dat', {'type', f_params{1}(1:end-4); 'parameter_set', 2; 'fitting_type', 'cmaes'});
              %vals = group_ml_results('adr-kymo-*_evol.dat', {'type', f_params{1}(1:end-4); 'parameter_set', 2; 'fitting_type', 'cmaes'});
              
              if (isempty(vals))
                warning('No adequat initial condition identified, skipping !');
                continue;
              end
              
              best = Inf;
              indx = 0;
              for j=1:size(vals{1,2}, 1)
                if (best > vals{1,2}{j,2}(end).score)
                  best = vals{1,2}{j,2}(end).score;
                  indx = j;
                end
              end

              tmp_p = vals{1,2}{indx,2}(end).params;
              tmp_p(1:4) = abs(tmp_p(1:4) .* vals{1,1}{1}.rescale_factor);
              tmp_p(5:end) = tmp_p(5:end) * vals{1,1}{1}.offset_scaling;
              s_params{2} = tmp_p;
            
            elseif (starts(s) == 'T')
              vals = group_ml_results('BestFits/adr-kymo-*_evol.dat', {'type', '1056-temps-all'; 'parameter_set', param_set(p); 'fitting_type', 'cmaes'});
              %vals = group_ml_results('adr-kymo-*_evol.dat', {'type', f_params{1}(1:end-4); 'parameter_set', 2; 'fitting_type', 'cmaes'});
              
              if (isempty(vals))
                warning('No adequat initial condition identified, skipping !');
                continue;
              end
              
              best = Inf;
              indx = 0;
              for j=1:size(vals{1,2}, 1)
                if (best > vals{1,2}{j,2}(end).score)
                  best = vals{1,2}{j,2}(end).score;
                  indx = j;
                end
              end

              tmp_p = vals{1,2}{indx,2}(end).params;
              s_params{2} = abs(tmp_p(1:4) .* vals{1,1}{1}.rescale_factor);
              %tmp_p(5:end) = tmp_p(5:end) * vals{1,1}{1}.offset_scaling;
              %s_params{2} = tmp_p;
            end

            all_params = [f_params; p_params; i_params; s_params; params];

            disp(' ');
            print_cmd(all_params{:}, [num2str(r) '/' num2str(repeats)]);
            disp([num2str(100*counts/ntotal) '%:' num2str(counts) '/' num2str(ntotal)]);
            disp(' ');
            counts = counts + 1;

            find_kymograph(all_params{:});
          end
        end
      end
    end
  end

  return;
end

function print_cmd(varargin)

  string = '';
  for i=1:length(varargin)
    string = [string num2str(varargin{i}) ','];
  end

  disp(string(1:end-1))

  return;
end
