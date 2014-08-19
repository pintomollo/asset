function perform_fitting(selection, roundit)

  if (nargin == 1)
    roundit = true;
  end

  all_files = {'good_24.txt'; 'good_20.txt'; 'good_13.txt'; 'good_c27d91.txt'; 'good_ani2.txt'};
  all_files = cellfun(@(x)(textread(x, '%s')), all_files, 'UniformOutput', false);

  params = {'config_fitting'; 'fit_kymo'; 'config_modeling'; 'custom_flow'; 'fixed_parameter'; []};
  init_noise = 0;
  fixed_parameter = [];
  user_param = [];

  switch_selection = selection;
  if (roundit)
    switch_selection = floor(selection);
  end

  % Default setting
  s_params = {};

  switch switch_selection
    case -4
      files = {'1056-24-all.mat'};
      repeats = 3;
      starts = 'G';
      param_set = 6;
      params{4} = 'goehring';
    case -3
      files = {'1056-24-all.mat'};
      repeats = 3;
      starts = 'G';
      param_set = 2;
      params{4} = 'goehring';
    case -2
      files = dir('1056-temps-all.mat');
      repeats = 3;
      init_noise = 0.01;
      starts = 'P';
      param_set = 15;
      %params{2} = 'refine_temp_indep';
      params{2} = 'fit_flows';
    case -1
      % All ani-2 and c27d9.1
      files = dir('1056-temps-all.mat');
      repeats = 1;
      init_noise = 0;
      starts = 'P';
      param_set = 15;
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
      starts = 'B';
      param_set = 2;
      params{2} = 'refine_fit';
    case 1
      % All 24°C
      files = dir('1056-24-*_.mat');
      repeats = 1;
      init_noise = 0.1;
      starts = 'A';
      param_set = 2;
    case 2
      % All 14°C
      files = dir('1056-14-*_.mat');
      repeats = 1;
      init_noise = 0.1;
      starts = 'A';
      param_set = 2;
    case 3
      % All 3°C, ani-2 and c27d9.1
      files = [dir('1056-3-*_.mat'); ...
               dir('1056-ani2-*_.mat'); ...
               dir('1056-c27d91-*_.mat')];
      repeats = 1;
      init_noise = 0.1;
      starts = 'A';
      param_set = 2;
    case 4
      perform_fitting(4.1,false);
      perform_fitting(4.2,false);
      perform_fitting(4.3,false);
      return;
    case 4.1
      files = {'1056-temps-all.mat'};
      repeats = 1;
      init_noise = [0 0.5 1];
      starts = 'T';
      param_set = 14;
      params{2} = 'refine_flow';
      params{4} = 'extended_model';
      fixed_parameter = [4];
    case 4.2
      files = {'1056-temps-all.mat'};
      repeats = 1;
      init_noise = [0 0.5 1];
      starts = 'T';
      param_set = [14 15];
      params{2} = 'refine_flow';
      params{4} = 'extended_model';
      fixed_parameter = [5];
    case 4.3
      files = {'1056-temps-all.mat'};
      repeats = 1;
      init_noise = [0 0.5 1];
      starts = 'T';
      param_set = 15;
      params{2} = 'refine_flow';
      params{4} = 'extended_model';
      fixed_parameter = [4 5];
    case 5
      files = {'1056-temps-all.mat'};
      repeats = 1;
      init_noise = [0 0.5 1];
      starts = 'T';
      param_set = [2 14 15];
      params{2} = 'refine_flow';
      params{4} = 'extended_model';
    case 6
      perform_fitting(6.1,false);
      perform_fitting(6.2,false);
      return;
    case 6.1
      files = {'1056-temps-all.mat'};
      repeats = 2;
      init_noise = [0 0.5 1];
      starts = 'O';
      param_set = [20 24];
      params{2} = 'fit_flows';
      params{4} = 'extended_model';
    case 6.2
      files = {'1056-temps-all.mat'};
      repeats = 2;
      init_noise = [0 0.5 1];
      starts = 'O';
      param_set = 24;
      params{2} = 'fit_temp_indep';
      params{4} = 'extended_model';
    case 6.3
      files = {'1056-temps-all.mat'};
      repeats = 2;
      init_noise = [0 0.5 1];
      starts = 'O';
      param_set = 15;
      params{2} = 'fit_flows';
      params{4} = 'extended_model';
    case 6.4
      files = {'1056-temps-all.mat'};
      repeats = 2;
      init_noise = [0 0.5 1];
      starts = 'O';
      param_set = 15;
      params{2} = 'fit_flows';
      params{4} = 'extended_model';
      fixed_parameter = [5];
    case 6.5
      files = {'1056-temps-all.mat'};
      repeats = 2;
      init_noise = [0 0.5 1];
      starts = 'O';
      param_set = 15;
      params{2} = 'fit_flows';
      params{4} = 'extended_model';
      fixed_parameter = [4 5];
    case 7
      perform_fitting(7.1,false);
      perform_fitting(7.2,false);
      return;
    case 7.1
      files = {'1056-temps-all.mat'};
      repeats = 2;
      init_noise = [0 0.5 1];
      starts = 'O';
      param_set = [20 24];
      params{2} = 'refine_flow';
      params{4} = 'extended_model';
    case 7.2
      files = {'1056-temps-all.mat'};
      repeats = 2;
      init_noise = [0 0.5 1];
      starts = 'O';
      param_set = 24;
      params{2} = 'refine_temp_indep';
      params{4} = 'extended_model';
    case 107
      files = {'simulation'};
      repeats = 10;
      init_noise = [0.1 0.25 0.5 1];
      starts = 'GA';
      param_set = 2;
    case 8
      perform_fitting(8.2,false);
      perform_fitting(8.3,false);
      perform_fitting(8.1,false);
      return;

    case 8.1
      files = dir('1056-all-all.mat');
      repeats = 1;
      init_noise = 0;
      %starts = 'O';
      starts = 'D';
      param_set = 30;
      params{2} = 'hessian';
      params{4} = 'final_model';
      fixed_parameter = [15:154];

      ps = load('all_offsets');
      s_params = {'init_pos'; ps.all_offsets(:,4).'};

    case 8.2
      files = dir('1056-all-all.mat');
      repeats = 1;
      init_noise = 0;
      %starts = 'O';
      starts = 'D';
      param_set = 30;
      params{2} = 'sensitivity';
      params{4} = 'final_model';
      fixed_parameter = [15:154];

      ps = load('all_offsets');
      s_params = {'init_pos'; ps.all_offsets(:,4).'};
    case 8.3
      files = dir('1056-all-all.mat');
      repeats = 1;
      init_noise = 0;
      %starts = 'O';
      starts = 1;
      param_set = 30;
      params{2} = 'refine_flow';
      params{4} = 'custom_flow_sample';
    case 9
      perform_fitting(9.1,false);
      perform_fitting(9.2,false);
      return;
    case 9.1
      %files = {'1056-all-all.mat'};
      files = [dir('1056-3-*_.mat'); ...
               dir('1056-14-*_.mat'); ...
               dir('1056-24-*_.mat'); ...
               dir('1056-ani2-*_.mat'); ...
               dir('1056-c27d91-*_.mat')];
      repeats = 1;
      init_noise = 0;
      starts = 'D';
      param_set = 0;
      params{2} = 'refine_fit';
    case 9.2
      %files = {'1056-all-all.mat'};
      files = [dir('1056-3-*_.mat'); ...
               dir('1056-14-*_.mat'); ...
               dir('1056-24-*_.mat'); ...
               dir('1056-ani2-*_.mat'); ...
               dir('1056-c27d91-*_.mat')];
      repeats = 1;
      init_noise = 0;
      starts = 'D';
      param_set = 0;
      params{2} = 'refine_fit';
      params{4} = 'goehring';
    case 9.3
      %files = {'1056-all-all.mat'};
      files = [dir('1056-3-*_.mat'); ...
               dir('1056-14-*_.mat'); ...
               dir('1056-24-*_.mat'); ...
               dir('1056-ani2-*_.mat'); ...
               dir('1056-c27d91-*_.mat')];
      repeats = 1;
      init_noise = 0;
      starts = 'D';
      param_set = 0;
      params{2} = 'refine_fit';
      params{4} = 'extended_model';
    case 9.4
      %files = {'1056-all-all.mat'};
      files = [dir('1056-3-*_.mat'); ...
               dir('1056-14-*_.mat'); ...
               dir('1056-24-*_.mat'); ...
               dir('1056-ani2-*_.mat'); ...
               dir('1056-c27d91-*_.mat')];
      repeats = 1;
      init_noise = 0;
      starts = 'D';
      param_set = 0;
      params{2} = 'refine_fit';
      params{4} = 'final_model';
    case 10
      files = cat(1, all_files{:});
      repeats = 1;
      init_noise = 0;
      starts = -2;
      param_set = 2;
      params{2} = 'refine_size';
      fixed_parameter = [1:4 6];

    case 11

      files = cat(1, all_files{:});
      %files = [dir('1056-3-*_.mat'); ...
      %         dir('1056-14-*_.mat'); ...
      %         dir('1056-24-*_.mat'); ...
      %         dir('1056-ani2-*_.mat'); ...
      %         dir('1056-c27d91-*_.mat')];
      repeats = 1;
      init_noise = 0;
      %starts = mod(selection, 10);
      starts = 'D';
      param_set = 0;
      params{2} = 'refine_fit';
      params{4} = 'extended_model';
      %fixed_parameter = [2];
    case 12
      files = {'1056-temps-all.mat'};
      repeats = 3;
      init_noise = 0;
      starts = 'T';
      param_set = 2;
      params{2} = 'fit_flows';
      params{4} = 'custom_model';
    case 13
      files = {'1056-temps-all.mat'};
      repeats = 3;
      init_noise = 0;
      starts = 'T';
      param_set = 15;
      params{2} = 'fit_flows';
      params{4} = 'custom_model';
    case 14
      files = {'1056-temps-all.mat'};
      repeats = 3;
      init_noise = 0;
      starts = 'T';
      param_set = 24;
      params{2} = 'fit_flows';
      params{4} = 'custom_model';
    case 15
      files = {'1056-temps-all.mat'};
      repeats = 3;
      init_noise = 0;
      starts = 'T';
      param_set = 24;
      params{2} = 'fit_temp_indep';
      params{4} = 'custom_model';

    case 16
      perform_fitting(16.1,false);
      perform_fitting(16.2,false);
      perform_fitting(16.3,false);
      return;
    case 16.1
      files = {'1056-temps-all.mat'};
      repeats = 1;
      init_noise = 0;
      starts = 3;
      param_set = 14;
      params{2} = 'refine_flow';
      params{4} = 'custom_model';
    case 16.2
      files = {'1056-temps-all.mat'};
      repeats = 1;
      init_noise = 0;
      starts = 4;
      param_set = 14;
      params{2} = 'refine_flow';
      params{4} = 'custom_model';
    case 16.3
      files = {'1056-temps-all.mat'};
      repeats = 1;
      init_noise = 0;
      starts = 8;
      param_set = 15;
      params{2} = 'refine_flow';
      params{4} = 'custom_model';
    case 17
      files = {'1056-temps-all.mat'};
      repeats = 3;
      init_noise = 0;
      starts = 3;
      param_set = 14;
      params{2} = 'fit_flows';
      params{4} = 'custom_model';
    case 18
      files = {'1056-temps-all.mat'};
      repeats = 3;
      init_noise = 0;
      starts = 4;
      param_set = 14;
      params{2} = 'fit_flows';
      params{4} = 'custom_model';
    case 19
      files = {'1056-temps-all.mat'};
      repeats = 3;
      init_noise = 0;
      starts = 8;
      param_set = 15;
      params{2} = 'fit_flows';
      params{4} = 'custom_model';
    case 20
      files = {'1056-temps-all.mat'};
      repeats = 1;
      init_noise = 0;
      starts = 'T';
      param_set = [14 20 14 20 14 20];
      params{2} = 'fit_flows';
      params{4} = 'custom_model';
    case 21
      files = {'1056-med-scale.mat'};
      repeats = 1;
      init_noise = 0;
      starts = 'P';
      param_set = 2;
      params{2} = 'fit_kymo';
    case 22
      files = {'1056-med-scale.mat'};
      repeats = 1;
      init_noise = 0;
      starts = 'D';
      param_set = 0;
      %params{2} = 'refine_border';
      params{2} = 'refine_size';
      params{4} = 'custom_model';
    case 23
      files = {'1056-med-scale.mat'};
      repeats = 5;
      init_noise = 0.05;
      starts = 'P';
      param_set = 0;
      params{2} = 'refine_size_flow';
      params{4} = 'custom_model';
      user_param = 10*(selection - switch_selection);
    case 24
      files = cat(1, all_files{:});
      %files = [dir('1056-3-*_.mat'); ...
      %         dir('1056-14-*_.mat'); ...
      %         dir('1056-24-*_.mat'); ...
      %         dir('1056-ani2-*_.mat'); ...
      %         dir('1056-c27d91-*_.mat')];
      repeats = 1;
      init_noise = 0;
      %starts = mod(selection, 10);
      starts = -6;
      param_set = 15;
      params{2} = 'refine_extended';
      params{4} = 'extended_model';
      fixed_parameter = [1:4 6:10];
    case 25
      files = {'1056-all-all.mat'};
      %files = [dir('1056-3-*_.mat'); ...
      %         dir('1056-14-*_.mat'); ...
      %         dir('1056-24-*_.mat'); ...
      %         dir('1056-ani2-*_.mat'); ...
      %         dir('1056-c27d91-*_.mat')];
      repeats = 1;
      init_noise = 0;
      %starts = mod(selection, 10);
      starts = 'D';
      param_set = 15;
      params{2} = 'refine_extended';
      params{4} = 'final_model';
      fixed_parameter = [2 4 145:149];
    case 26
      files = {'1056-all-all.mat'};
      %files = [dir('1056-3-*_.mat'); ...
      %         dir('1056-14-*_.mat'); ...
      %         dir('1056-24-*_.mat'); ...
      %         dir('1056-ani2-*_.mat'); ...
      %         dir('1056-c27d91-*_.mat')];
      repeats = 1;
      init_noise = 0;
      %starts = mod(selection, 10);
      starts = 'D';
      param_set = 15;
      params{2} = 'refine_extended';
      params{4} = 'final_model';
      fixed_parameter = [2 4 149];
    case 27
      files = {'1056-all-all.mat'};
      %files = [dir('1056-3-*_.mat'); ...
      %         dir('1056-14-*_.mat'); ...
      %         dir('1056-24-*_.mat'); ...
      %         dir('1056-ani2-*_.mat'); ...
      %         dir('1056-c27d91-*_.mat')];
      repeats = 1;
      init_noise = 0;
      %starts = mod(selection, 10);
      starts = 'D';
      param_set = 15;
      params{2} = 'fit_extended';
      params{4} = 'final_model';
      fixed_parameter = [2 4 145:149];
    case 28
      files = {'1056-all-all.mat'};
      %files = [dir('1056-3-*_.mat'); ...
      %         dir('1056-14-*_.mat'); ...
      %         dir('1056-24-*_.mat'); ...
      %         dir('1056-ani2-*_.mat'); ...
      %         dir('1056-c27d91-*_.mat')];
      repeats = 1;
      init_noise = 0;
      %starts = mod(selection, 10);
      starts = 'D';
      param_set = 15;
      params{2} = 'refine_extended';
      params{4} = 'final_model';
      fixed_parameter = [2 4 5:147];
    case 29
      files = {'1056-all-all.mat'};
      %files = [dir('1056-3-*_.mat'); ...
      %         dir('1056-14-*_.mat'); ...
      %         dir('1056-24-*_.mat'); ...
      %         dir('1056-ani2-*_.mat'); ...
      %         dir('1056-c27d91-*_.mat')];
      repeats = 1;
      init_noise = 0;
      %starts = mod(selection, 10);
      starts = 'D';
      param_set = 15;
      params{2} = 'refine_extended';
      params{4} = 'final_model';
      fixed_parameter = [2 4 5:144 149];
    case 30
      files = {'1056-24-all.mat'};
      repeats = 1;
      init_noise = 0;
      starts = [10:13];
      param_set = 2;
      params{2} = 'refine_fit';
      params{4} = 'custom_flow';
    case 31
      files = {'1056-med-all.mat'};
      repeats = 1;
      init_noise = 0;
      starts = 'D';
      param_set = 2;
      params{2} = 'refine_fit';
      params{4} = 'final_model';
      fixed_parameter = [2 4];
    case 32
      files = {'1056-med-all.mat'};
      repeats = 1;
      init_noise = 0;
      starts = 'D';
      param_set = 2;
      %params{2} = 'refine_border';
      params{2} = 'refine_size';
      params{4} = 'final_model';
      fixed_parameter = [2 4];
    case 33
      files = {'1056-med-all.mat'};
      repeats = 1;
      init_noise = 0;
      starts = 'D';
      param_set = 2;
      params{2} = 'refine_fit';
      params{4} = 'final_model';
    case 34
      files = {'1056-med-all.mat'};
      repeats = 1;
      init_noise = 0;
      starts = 'D';
      param_set = 2;
      %params{2} = 'refine_border';
      params{2} = 'refine_size';
      params{4} = 'final_model';
    case 35
      files = {'1056-med-all.mat'};
      repeats = 2;
      init_noise = 0;
      starts = 'D';
      param_set = 2;
      params{2} = 'fit_kymo';
      params{4} = 'final_model';
      fixed_parameter = [2 4];
    case 36
      files = {'1056-med-all.mat'};
      repeats = 2;
      init_noise = 0;
      starts = 'D';
      param_set = 2;
      %params{2} = 'refine_border';
      params{2} = 'fit_size';
      params{4} = 'final_model';
      fixed_parameter = [2 4];

    case 40
      files = {'1056-size-all.mat'};
      repeats = 5;
      init_noise = 0.05;
      starts = 'D';
      param_set = 2;
      params{2} = 'refine_fit';
      params{4} = 'extended_model';
    case 41
      files = {'1056-med-scale.mat'};
      repeats = 3;
      init_noise = 0.05;
      starts = 'P';
      param_set = 2;
      params{2} = 'fit_kymo';
      params{4} = 'extended_model';
      fixed_parameter = [2 4];
    case 42
      files = {'1056-med-scale.mat'};
      repeats = 5;
      init_noise = 0.05;
      starts = 'P';
      param_set = 2;
      params{2} = 'refine_fit';
      params{4} = 'extended_model';
      %fixed_parameter = [2 4];
    case 43
      perform_fitting(43.1,false);
      perform_fitting(43.2,false);
      return;
    case 43.1
      files = {'1056-temps-all.mat'};
      repeats = 3;
      init_noise = 0;
      starts = 'D';
      param_set = 15;
      params{2} = 'fit_kymo';
      params{4} = 'full_model';
      fixed_parameter = [2 4];
    case 43.2
      files = {'1056-temps-all.mat'};
      repeats = 3;
      init_noise = 0;
      starts = 'D';
      param_set = 15;
      params{2} = 'refine_fit';
      params{4} = 'full_model';
    case 44
      perform_fitting(44.1,false);
      perform_fitting(44.2,false);
      return;
    case 44.1
      files = {'1056-temps-all.mat'};
      repeats = 3;
      init_noise = 0;
      starts = 'D';
      param_set = 15;
      params{2} = 'refine_fit';
      params{4} = 'full_model';
      fixed_parameter = [2 4];
    case 44.2
      files = {'1056-temps-all.mat'};
      repeats = 3;
      init_noise = 0;
      starts = 'D';
      param_set = 15;
      params{2} = 'fit_kymo';
      params{4} = 'full_model';
    case 45
      perform_fitting(45.1,false);
      perform_fitting(45.2,false);
      return;
    case 45.1
      files = {'1056-temps-all.mat'};
      repeats = 3;
      init_noise = 0;
      starts = 'D';
      param_set = 15;
      params{2} = 'fit_kymo';
      params{4} = 'final_model';
      fixed_parameter = [2 4];
    case 45.2
      files = {'1056-temps-all.mat'};
      repeats = 3;
      init_noise = 0;
      starts = 'D';
      param_set = 15;
      params{2} = 'refine_fit';
      params{4} = 'final_model';
      fixed_parameter = [2 4 9];
    case 46
      perform_fitting(46.1,false);
      perform_fitting(46.2,false);
      return;
    case 46.1
      files = {'1056-temps-all.mat'};
      repeats = 3;
      init_noise = 0;
      starts = 'D';
      param_set = 15;
      params{2} = 'refine_fit';
      params{4} = 'final_model';
      fixed_parameter = [2 4 ];
    case 46.2
      files = {'1056-temps-all.mat'};
      repeats = 3;
      init_noise = 0;
      starts = 'D';
      param_set = 15;
      params{2} = 'fit_kymo';
      params{4} = 'final_model';
    case 47
      files = {'1056-temps-all.mat'};
      repeats = 3;
      init_noise = 0.01;
      starts = 'D';
      param_set = 15;
      params{2} = 'refine_flow';
      params{4} = 'full_model';
      fixed_parameter = [2 4];
    case 48
      files = {'1056-all-all.mat'};
      repeats = 3;
      init_noise = 0.01;
      starts = 'D';
      param_set = 15;
      params{2} = 'refine_fit';
      params{4} = 'full_model';
      fixed_parameter = [2 4];
    case 49
      files = cat(1, all_files{:});
      repeats = 1;
      init_noise = 0;
      starts = 'D';
      param_set = 0;
      params{2} = 'refine_fit';
      params{4} = 'extended_model';

    otherwise
      warning('Choose a fitting group between 1 and 7');

      return;
  end

  params{6, 1} = fixed_parameter;

  ntotal = length(files)*repeats*length(init_noise)*length(starts)*length(param_set);
  counts = 1;

  for s = 1:length(starts)
    switch starts(s)
      case 'A'
      % Average of all bests
        s_params = {'init_pos'; [0.00769 2.197 0.0314 2.202]};
      case 'P'
      % Provided

        %s_params = {'init_pos'; [0.00154 2.2569 0.0078 2.0203 -2.9900 9.6900 0 0.1599 0 0.6277 1.6908 0.8265]};
        % refine_flow for 25
        %s_params = {'init_pos'; [0.00154 2.2569 0.0078 2.0203 -2.9900 9.6900 0 0.1599 0.1599 0.1599 0 0.6277 1.6908 0.8265]};
        %s_params = {'init_pos'; [0.00769 2.197 0.0314 2.202 -2.9900 9.6900 0 0.1599 0 0.83 0.8265]};

        %s_params = {'init_pos'; [0.0116 2.1571 0.0658 2.1871 21 2.2400 28.0400]};
        %if (user_param == 0)
        %  s_params = {'init_pos'; [18.766 -1.9869 33.482 0.856 1.5438]};
        %else
        %  s_params = {'init_pos'; [18.766 -1.9869 33.482 0.856 user_param]};
        %end

        %% Revisions....
        s_params = {'init_pos'; [0.002519 2.2033 0.015445 2.2487 12.36 -15 30.35]};
        %s_params = {'init_pos'; [0.002721 2.1571 0.016588 2.1871 12.73 -6.99 38.82]};

        %refine_temp_indep
        %s_params = {'init_pos'; [0.00154 2.2569 0.0078 2.0203 -2.9900 9.6900 0 0.8566 0.8265 1.0889 0.8265 0.6127 1.7139 0.8265]};

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
          %ps = load('full_params');
          %s_params = {'init_pos'; ps.params(starts(s),:)};
          ps = load('all_offsets');
          tmp_data = load('temp_params.mat');

          switch abs(starts(s))
            case 1
              s_params = {'init_pos'; [ps.all_offsets.' 1.441]};
            case 2
              s_params = {'init_pos'; [0.0116 2.1571 0.0658 2.1871 ps.all_offsets.' 1.441]};
            case {3, 4}
              good_indx = (tmp_data.param_set(:,1) == param_set);
              s_params = {'init_pos'; [0.0116 2.1571 0.0658 2.1871 tmp_data.params{good_indx}(8:end) 1.441]};
              params{6, 1} = [2 4 (length(s_params{2})+3)];

              if (starts == 3)
                s_params{2}(end-3) = 0;
                params{6, 1} = [params{6,1}(1:2) (params{6,1}(3)-[3 0])];
              else
                s_params{2}(end-2) = 0;
                params{6, 1} = [params{6,1}(1:2) (params{6,1}(3)-[2 0])];
              end
            case 5
              good_indx = (tmp_data.param_set(:,1) == param_set);
              s_params = {'init_pos'; [0.0116 2.1571 0.0658 2.1871 tmp_data.params{good_indx}(8:end) 1.441]};
              params{6, 1} = [2 4 (length(s_params{2})+3)];

              %s_params{2}([end-4 end-3]) = 0;
              %params{6, 1} = [params{6,1}(1:2) (params{6,1}(3)-[4 3 0])];
            case {6, 7}
              s_params = {'init_pos'; [0.002587 2.1571 0.01482 2.1871 ps.all_offsets.' 0.3566 0 0.9734 0.8755 1.441]};
            case 8
              good_indx = (tmp_data.param_set(:,1) == param_set);
              s_params = {'init_pos'; [0.0116 2.1571 0.0658 2.1871 tmp_data.params{good_indx}(8:end) 1.441]};
              params{6, 1} = [2 4 (length(s_params{2})+[0 3])];
              s_params{2}(end-3) = 0;

            case 10
              s_params = {'init_pos'; [0.00769 2.197 0.0314 2.202]};
              params{6, 1} = [1:4 6:length(s_params{2})+1];
              params{2} = 'refine_fit';
            case 11
              s_params = {'init_pos'; [0.0116 2.1571 0.0658 2.1871]};
              params{6, 1} = [1:4 6:length(s_params{2})+1];
              params{2} = 'refine_fit';
            case 12
              s_params = {'init_pos'; [0.0116 2.1571 0.0658 2.1871 1.4575]};
              params{6, 1} = [1:4 6:length(s_params{2})+1];
              params{2} = 'refine_extended';
            case 13
              s_params = {'init_pos'; [0.002535 2.1571 0.014496 2.1871 0.2523 0.0739 0.9062 0.9409 1.4575]};
              params{6, 1} = [1:4 6:length(s_params{2})+1];
              params{2} = 'refine_flow';
              param_set = 15;
            otherwise
              s_params = {'init_pos'; ps.all_offsets.'};
          end

          orig_s = s_params;
        %else
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

              tmp_opts = load_parameters('modeling', params{4});
              tmp_vals = tmp_opts.reaction_params([3:4],:);
              tmp_vals = [tmp_vals(:).' -3.9800 9.4600 0.3900];
              fixed_parameter = fixed_parameter + 4;

              switch param_set(p)
                case 15
                  tmp_vals = [tmp_vals 0.2523 0.0739 0.9062 0.9409 1.4575];
                  tmp_vals(fixed_parameter) = 0;
                case 20
                  tmp_vals = [tmp_vals 0.4059 1.6429 0.9409 1.4575];
                case 24
                  if (params{2}(end) == 'p')
                    tmp_vals = [tmp_vals 0.78832 0.9127 1.1439 1.0197 0.4059 1.6429 0.940 1.4575];
                  else
                    tmp_vals = [tmp_vals 0.2523 0.0739 0.4059 1.6429 0.9409 1.4575];
                  end
              end
              s_params{2} = tmp_vals;
              params{6, 1} = [2 4 fixed_parameter length(tmp_vals)];

            elseif (starts(s) == 'T')

              tmp_data = load('temp_params.mat');
              tmp_opts = load_parameters('modeling', params{4});

              good_indx = (tmp_data.param_set(:,1) == param_set(p));
              if (any(good_indx))
                if (sum(good_indx) > 1)
                  ind = fixed_parameter;
                  if (isempty(ind))
                    ind = 0;
                  elseif (numel(ind) > 1)
                    ind = 10*ind(1) + ind(2);
                  end

                  good_indx = good_indx & (tmp_data.param_set(:,3) ~= (params{2}(end) == 'p') & (tmp_data.param_set(:,2) == ind));
                end

                tmp_vals = tmp_opts.reaction_params([3:4],:);
                s_params{2} = [tmp_vals(:).' tmp_data.params{good_indx, 1}(5:end)];
                params{6, 1} = [2 4 fixed_parameter+4];
              else
                warning('No adequat initial condition identified, skipping !');
                continue;
              end
            elseif (starts(s) < 0)
              switch abs(starts)
                case 1
                  s_params{2,1} = orig_s{2,1}([f end]);
                case 2
                  s_params{2,1} = orig_s{2,1}([1:4 (f+4) end]);
                case 6
                  s_params{2,1} = orig_s{2,1}([1:4 (f+4) end-4:end]);
              end
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
