function perform_fitting(selection)

  params = {'config_fitting'; 'fit_kymo'; 'config_modeling'; 'custom_flow'};
  init_noise = 0;
  
  switch selection
    case -1
      % All ani-2 and c27d9.1
      files = dir('1056-temps-all.mat');
      repeats = 1;
      init_noise = 0;
      starts = 'P';
      param_set = 24;
      %params{2} = 'refine_temp_indep';
      params{2} = 'refine_flow';
      %params{2} = 'refine_fit';
    case 0
      % All ani-2 and c27d9.1
      files = [dir('1056-3-*_.mat'); ...
               dir('1056-14-*_.mat'); ...
               dir('1056-24-*_.mat'); ...
               dir('1056-ani2-*_.mat'); ...
               dir('1056-c27d91-*_.mat')];
      repeats = 1;
      init_noise = 0;
      starts = 'A';
      param_set = 2;
      params{2} = 'refine_fit';
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
      repeats = 1;
      init_noise = 0;
      starts = 'B';
      param_set = 2;
      params{2} = 'refine_fit';
    case 5
      files = {'1056-temps-all.mat'};
      repeats = 2;
      init_noise = 0.1;
      starts = 'AB';
      param_set = [12 13 23];
      params{2} = 'refine_fit';
    case 6
      files = {'1056-temps-all.mat'};
      repeats = 2;
      init_noise = 0.1;
      starts = 'AB';
      param_set = [24 25];
      params{2} = 'refine_fit';
    case 7
      files = {'simulation'};
      repeats = 10;
      init_noise = [0.1 0.25 0.5 1];
      starts = 'GA';
      param_set = 2;
    case 8
      perform_fitting(8.2);
      perform_fitting(8.3);
      perform_fitting(8.1);
      return;

    case 8.1
      files = dir('1056-temps-all.mat');
      repeats = 1;
      init_noise = 0;
      starts = 'O';
      param_set = 30;
      params{2} = 'hessian';
      params{4} = 'custom_flow_sample';
    case 8.2
      files = dir('1056-temps-all.mat');
      repeats = 1;
      init_noise = 0;
      starts = 'O';
      param_set = 30;
      params{2} = 'sensitivity';
      params{4} = 'custom_flow_sample';
    case 8.3
      files = dir('1056-temps-all.mat');
      repeats = 1;
      init_noise = 0;
      starts = 'O';
      param_set = 30;
      params{2} = 'refine_flow';
      params{4} = 'custom_flow_sample';
    case 9
      files = {'1056-temps-all.mat'};
      repeats = 1;
      init_noise = 0;
      starts = 'B';
      param_set = [2 12 13 20 23 24 25];
      params{2} = 'refine_fit';
    case 10
      files = {'1056-temps-all.mat'};
      repeats = 1;
      init_noise = 0;
      starts = 'B';
      param_set = [2 12 13 20 23 24 25];
      params{2} = 'refine_flow';
    case 11
      files = {'1056-all-all.mat'};
      repeats = 1;
      init_noise = 0;
      starts = mod(selection, 10);
      param_set = 24;
      params{2} = 'refine_flow';
    case 12
      files = {'1056-all-all.mat'};
      repeats = 1;
      init_noise = 0;
      starts = mod(selection, 10);
      param_set = 24;
      params{2} = 'refine_flow';
    case 13
      files = {'1056-all-all.mat'};
      repeats = 1;
      init_noise = 0;
      starts = mod(selection, 10);
      param_set = 24;
      params{2} = 'refine_flow';
    case 14
      files = {'1056-temps-all.mat'};
      repeats = 2;
      init_noise = 0;
      starts = 'T';
      param_set = 24;
      params{2} = 'fit_flows';
    case 15
      files = {'1056-temps-all.mat'};
      repeats = 2;
      init_noise = 0.1;
      starts = 'T';
      param_set = 24;
      params{2} = 'fit_temp_indep';
    case 16
      files = {'1056-temps-all.mat'};
      repeats = 2;
      starts = 'AB';
      param_set = 24;
      params{2} = 'refine_temp_indep';

    otherwise
      warning('Choose a fitting group between 1 and 7');

      return;
  end

  ntotal = length(files)*repeats*length(init_noise)*length(starts)*length(param_set);
  counts = 1;

  for s = 1:length(starts)
    switch starts(s)
      case 'A'
      % Average of all bests
        s_params = {'init_pos'; [0.00286 2.22 0.0153 2.31]};
      case 'P'
      % Provided

        s_params = {'init_pos', [0.00154 2.2569 0.0078 2.0203 -2.9900 9.6900 0 0.1599 0 0.6277 1.6908 0.8265]};
        % refine_flow for 25
        %s_params = {'init_pos', [0.00154 2.2569 0.0078 2.0203 -2.9900 9.6900 0 0.1599 0.1599 0.1599 0 0.6277 1.6908 0.8265]};
        %refine_temp_indep
        %s_params = {'init_pos', [0.00154 2.2569 0.0078 2.0203 -2.9900 9.6900 0 0.8566 0.8265 1.0889 0.8265 0.6127 1.7139 0.8265]};

       %s_params = {'init_pos'; [0.00134 1.7874 0.0091 2.4993 -5.9900 4.4900 0 0.6652 0.9319 1.2515 1.0152 0.3760 2.2994 0.9848]};
       %s_params = {'init_pos'; [0.00083 2.3309 0.0050 2.4398 -5.9900 8.1600 0 0 0 1.0161 1.6961 0.8069]};
       %s_params = {'init_pos'; [0.00171 1.8601 0.0104 2.2477 -5.0000 8.8100 -2.0000 0.0195 0.0195 0.8283 1.7989]};
      case 'B'
      % Best value found
        s_params = {'start_with_best'; true};
      case {'O', 'T'}
      % Current best found during optimization in BestFits
        s_params = {'init_pos'; []};
      otherwise
        if (isnumeric(starts))
          ps = load('full_params');
          s_params = {'init_pos'; ps.params(starts(s),:)};
        else
        % Default setting
          s_params = {};
        end
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
              %{
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
              %}
              s_params{2} = [0.28 0.00857 0.0054 0.00154 2.2569 1.56 0.17353 67.5 0.15 0.0472 0.0073 0.0078 2.0203 1 -2.9900 9.6900 0 0.1599 0 0.6277 1.6908 0.8265];

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
