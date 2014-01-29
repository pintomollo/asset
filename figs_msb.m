function figs_msb(num)

  if (nargin == 0)
    num = 1;
  end
  opts_expansion = load_parameters(get_struct('ASSET'), 'domain_expansion.txt');
  colors = 'gbrcmyk';
  black = ones(1,3)*37/255;

  thresh = [0 0.775];
  nframes = 25;
  ntails = 20;

  switch num
    case 0

      all_files = {'good_24.txt'; 'good_20.txt'; 'good_13.txt'; 'good_c27d91.txt'; 'good_ani2.txt'};
      all_data = cell(5, 3);
      all_data(:,1) = all_files(:);

      [time, names] = get_manual_timing();

      for f = 1:length(all_files)
        files = textread(all_files{f}, '%s');
        nfiles = length(files);

        all_params = NaN(nfiles, 8);
        all_profiles = cell(nfiles, 3);
        for i=1:nfiles
          load(files{i});
          good_time = ismember(names, files{i});

          [profile, center, max_width, cell_width, path] = get_profile(mymovie, nframes, opts);
          all_params(i,:) = [2*max_width, cell_width, mymovie.metadata.axes_length_3d(:).' time(good_time, :)];
          all_profiles{i,1} = profile;
          all_profiles{i,2} = path;
          all_profiles{i,3} = files{i};

          disp([num2str(i) '/' num2str(nfiles)]);
        end
        all_data{f,2} = all_params;
        all_data{f,3} = all_profiles;
      end

      temperatures = [24 20 13 23 23];

      save('data_expansion', 'all_data', 'temperatures');

      keyboard

    case 0.1
      vals = group_ml_results('BestFits/adr-kymo-*_evol.dat', {'type'}, {'parameter_set', 2; 'fit_flow', false; 'extrapol_z', true; 'rescale_length_only', true; 'scale_each_egg', true});
      vals = extract_model_parameters(vals, true);
      nfits = size(vals,1);

      all_files = {'good_24.txt'; 'good_20.txt'; 'good_13.txt'; 'good_c27d91.txt'; 'good_ani2.txt'};
      all_data = cell(size(all_files, 1)+1, 3);

      for i=1:length(all_files)
        all_data{i,1} = all_files{i};
        all_data{i,3} = textread(all_files{i}, '%s');
        all_data{i,2} = NaN(size(all_data{i,3},1), 9);
      end
      all_data{end,1} = 'averages';
      all_data{end,2} = NaN(0,9);

      ngroups = size(all_data, 1);

      data = load('data_expansion');

      for i=1:nfits
        group_indx = -1;
        sub_indx = -1;
        curr_params = vals{i,1}{1};
        curr_name = curr_params.type;
        curr_vals = vals{i,2};

        if (strncmp(curr_params.fitting_type, 'sample', 6))
          continue;
        end

        for j=1:ngroups
          goods = ismember(all_data{j,3}, [curr_name '.mat']);
          if (any(goods))
            group_indx = j;
            sub_indx = find(goods, 1);

            break;
          end
        end

        if (group_indx < 0 && strncmp('all', curr_name(end-2:end), 3))
          group_indx = ngroups;
          all_data{group_indx, 2}(end+1,:) = NaN;
          sub_indx = size(all_data{group_indx, 2}, 1);

          all_data{group_indx, 3}{sub_indx} = [curr_name '.mat'];
        end

        if (group_indx > 0)
          curr_score = Inf;
          curr_pos = [];
          for j=1:size(curr_vals, 1)
            if (curr_vals{j,2}.score < curr_score)
              curr_score = curr_vals{j,2}.score;
              %curr_pos = abs(curr_vals{j,2}(end).params(1:4)) .* curr_params.rescale_factor;
              %curr_pos = [curr_pos curr_vals{j,2}(end).params(5)];
              curr_pos = curr_vals{j,2}.params.rate;
            end
          end

          if (isfinite(curr_score))
            if (group_indx ~= ngroups)
              good_group = find(ismember(data.all_data(:,1), all_data{group_indx,1}));
              good_sub = find(ismember(data.all_data{good_group,3}(:,end), all_data{group_indx,3}{sub_indx}));

              egg_size = data.all_data{good_group, 2}(good_sub, 3:5);
              npts = length(data.all_data{good_group, 3}{good_sub, 2});
            else
              egg_size = NaN(1, 3);
              movie = load(curr_name);
              if (iscell(movie.ground_truth))
                npts = 0;
                for g=1:length(movie.ground_truth)
                  npts = npts + size(movie.ground_truth{g}, 1);
                end
              else
                npts = size(movie.ground_truth, 1);
              end
            end

            prev_score = all_data{group_indx,2}(sub_indx, 1);
            if (isnan(prev_score) || curr_score < prev_score)
              all_data{group_indx,2}(sub_indx, :) = [curr_score, curr_pos, egg_size, npts];
            end
          end
        end
        disp(curr_name);
      end

      save('data_fitting', 'all_data');
      keyboard;

    case 0.2
      load('data_expansion')

      opts = get_struct('modeling');
      opts = load_parameters(opts, 'goehring.txt');
      opts = load_parameters(opts, 'custom_flow.txt');

      flow = opts.advection_params;
      if (size(flow, 1) ~= opts.nparticles)
        [X, Y] = meshgrid([1:size(flow, 2)], 1+([0:opts.nparticles-1]*(size(flow, 1)-1)/(opts.nparticles-1)).');
        flow = bilinear_mex(flow, X, Y, [2 2]);
      end

      for f = 1:size(all_data, 1)
        for i=1:size(all_data{f,2}, 1)
          egg_size = all_data{f,2}(i, 3:5).';
          cell_width = all_data{f,2}(i,2);

          egg_size(1) = ellipse_circum(egg_size, cell_width, true);

          opts.axes_length = egg_size;

          opts.reaction_params(end-1,:) = surface2volume(opts.axes_length);
          opts.reaction_params(end, :) = 0.5*ellipse_circum(opts.axes_length);

          opts.boundaries = [0 opts.reaction_params(end,1)];
          opts.x_step = diff(opts.boundaries)/(opts.nparticles-1);

          x0 = opts.init_func(opts, false);

          ml_params = [opts.diffusion_params; ...
                        opts.reaction_params];

          [domain, t_pos] = simulate_model_mix(x0, ml_params, opts.x_step, opts.tmax*0.75, opts.time_step, opts.output_rate, flow, opts.user_data, opts.max_iter);

          domain = domain((end/2)+1:end, :).';
          domain = [domain domain(:,end-1:-1:1)];

          [profile, center, max_width, cell_width, path] = get_profile(domain, nframes);
          norig = length(all_data{f,3}{i,1})-1;
          nprofile = length(profile)-1;
          profile = interp1([0:nprofile], profile, [0:norig]*nprofile/norig);

          all_data{f,2}(i,1) = 2*max_width * opts.x_step/opts_expansion.quantification.resolution;
          all_data{f,3}{i,1} = profile;
          all_data{f,3}{i,2} = path;

          disp([num2str(i) '/' num2str(size(all_data{f,2},1))]);
        end
      end

      all_simul = all_data;

      save('simul_expansion', 'all_simul', 'temperatures');

      keyboard

    case 0.3
      load('data_expansion')

      temps = unique(temperatures);
      ntemps = length(temps);
      good_visc = (temps~=20);
      nvisc = sum(good_visc);

      kB = 8.6173324e-5;
      C2K = 273.15;

      load('full_params');
      params = params(:,[1:4 end-5:end]);

      %new_vals = [1.55:0.05:1.95].';
      %params(1,3) = 0.01;

%      params = repmat(params(1,:), length(new_vals), 1);
%      params = [params new_vals];

%      params = [interp1([1 5], params([1 3], 1:4), [1:5]) repmat(params(1,5:end), 5, 1)];
%      params = params(3:end,:);
%      params = params(4,:);
%      new_vals = [0.8:0.1:1.2].';
%      params = repmat(params, length(new_vals), 1);
%      params(:,3) = params(:,3) .* new_vals;

      params = params(1:3,:);

      params(1,1:4) = [0.19 1 2 2];
      params(2,1:4) = [0.00556 2.161 0.0307 2.013];
      params(3,1:4) = [0.00596 2.038 0.0304 1.835];

      %params(end) = 0.8;
      %params(end-4) = 0.05;

      params(:, end-5:end) = repmat([0 0 1 1 1 1], size(params,1), 1);

      full_params = {};
      for i=1:size(params, 1)
        full_params{end+1} = params(i,:);
      end

      opts = get_struct('modeling');
      opts = load_parameters(opts, 'goehring.txt');
      opts = load_parameters(opts, 'custom_flow.txt');

      flow = opts.advection_params;
      if (size(flow, 1) ~= opts.nparticles)
        [X, Y] = meshgrid([1:size(flow, 2)], 1+([0:opts.nparticles-1]*(size(flow, 1)-1)/(opts.nparticles-1)).');
        flow = bilinear_mex(flow, X, Y, [2 2]);
      end

      orig_opts = opts;

      for p = 1:length(full_params)

        tmp_opts = orig_opts;
        params = full_params{p};

        switch length(params)
          % Parameter set 24
          case 10
            tmp_opts.reaction_params([3 4 10 11]) = params(1:4);
            more_params = params(5:end);

          % Parameter set 24 with test of rho
          case 11
            tmp_opts.reaction_params([3 4 10 11]) = params(1:4);
            tmp_opts.reaction_params(5) = params(end);
            more_params = params(5:end-1);

          % Parameter set 30
          case 20

            tmp = [tmp_opts.diffusion_params; tmp_opts.reaction_params];
            tmp(1:14) = params(1:14);

            tmp_opts.diffusion_params = tmp(1,:);
            tmp_opts.reaction_params = tmp(2:end,:);

            more_params = params(15:end);

          otherwise
            warning('Unknown parameter length.')
            continue;
        end

        curr_flow_scale = abs(more_params(end));
        more_params = more_params(1:end-1);

        curr_visc = ones(1,ntemps);
        curr_visc(good_visc) = abs(more_params(end-nvisc+1:end));
        more_params = more_params(1:end-nvisc);

        E = ones(3,2)*abs(more_params(end-1));
        flow_E = abs(more_params(end));

        for f = 1:size(all_data, 1)

          diff_ratio = curr_visc(temps==temperatures(f))*((temperatures(f)+C2K) / (opts.reaction_temperature+C2K));
          ratio = exp(-(E/kB)*((1/(temperatures(f)+C2K)) - (1/(opts.reaction_temperature+C2K))));
          flow_ratio = exp(-(flow_E/kB)*((1/(temperatures(f)+C2K)) - (1/(opts.flow_temperature+C2K))));

          temp_scaling = [ones(1,2)*diff_ratio; ratio; ones(4,2)];
          temp_flow_scale = curr_flow_scale * flow_ratio;

          for i = 1:size(all_data{f,2}, 1)

            opts = tmp_opts;

            egg_size = all_data{f,2}(i, 3:5).';
            cell_width = all_data{f,2}(i,2);

            egg_size(1) = ellipse_circum(egg_size, cell_width, true);

            opts.axes_length = egg_size;

            opts.reaction_params(end-1,:) = surface2volume(opts.axes_length);
            opts.reaction_params(end, :) = 0.5*ellipse_circum(opts.axes_length);

            opts.boundaries = [0 opts.reaction_params(end,1)];
            opts.x_step = diff(opts.boundaries)/(opts.nparticles-1);

            ml_params = [opts.diffusion_params; ...
                          opts.reaction_params];
            ml_params = ml_params .* temp_scaling;

            opts.diffusion_params = ml_params(1, :);
            opts.reaction_params = ml_params(2:end, :);

            x0 = opts.init_func(opts, true, false);

            [domain, t_pos] = simulate_model_real(x0, ml_params, opts.x_step, opts.tmax*0.75, opts.time_step, opts.output_rate, flow*temp_flow_scale, opts.user_data, opts.max_iter);

            domain = domain((end/2)+1:end, :).';
            domain = [domain domain(:,end-1:-1:1)];

            [profile, center, max_width, cell_width, path] = get_profile(domain, nframes);

            %{
            figure;subplot(1,2,1);
            imagesc(domain);
            subplot(1,2,2);
            plot(path);
            mtit(num2str(i))
            (find(path > 0, 1) - find(path > 0.775, 1))
            keyboard
            %}

            norig = length(all_data{f,3}{i,1})-1;
            nprofile = length(profile)-1;
            profile = interp1([0:nprofile], profile, [0:norig]*nprofile/norig);

            all_data{f,2}(i,1) = min(2*max_width * opts.x_step/opts_expansion.quantification.resolution, all_data{f,2}(i,2));
            all_data{f,3}{i,1} = profile;
            all_data{f,3}{i,2} = path;

            disp([num2str(i) '/' num2str(size(all_data{f,2},1)) ':' num2str(p)]);
          end
          %keyboard
        end

        all_simul{p} = all_data;

      end

      save('simul_optimized', 'all_simul', 'temperatures', 'full_params');
      figs_msb(8.1);
      figs_msb(8.2);

      keyboard
    case 1
      files = dir('stainings/*.txt');
      nfiles = length(files);
      names = cell(nfiles, 1);

      objective = 'Mean';
      %objective = 'Max';

      all_vals = NaN(0,4);
      for i=1:nfiles
        tmp = importdata(['stainings/' files(i).name]);
        %keyboard
        %good_vals = ismember(tmp.colheaders, 'Mean') | ismember(tmp.colheaders, 'Median');

               %tmp.data(:, ismember(tmp.colheaders, 'Max')), ...
               %tmp.data(:, ismember(tmp.colheaders, 'Mode')), ...
               %2*tmp.data(:, ismember(tmp.colheaders, 'StdDev')), ...
        tmp = [tmp.data(:, ismember(tmp.colheaders, 'Slice')), ...
               tmp.data(:, ismember(tmp.colheaders, objective)), ...
               tmp.data(:, ismember(tmp.colheaders, 'BX')), ...
               tmp.data(:, ismember(tmp.colheaders, 'Width'))];

        %tmp = tmp.data(:, good_vals);
        %tmp = [tmp(1:2:end,:) - tmp(2:2:end,:) tmp(1:2:end,:) tmp(2:2:end,:)];

        slices = unique(tmp(:,1)).';
        all_pts = NaN(length(slices), 3);

        for j=slices
          pts = tmp(tmp(:,1)==j, 2:end);
          [junk, ant] = min(pts(:,2));
          [junk, post] = max(pts(:,2)+pts(:,3));

          bkg = [1:size(pts,1)];
          bkg = bkg(bkg ~= ant & bkg ~= post);
          bkg = sort(pts(bkg, 1));

          dint = diff(bkg);

          if (pts(ant,1) > bkg(1))
            all_pts(j,:) = [(pts(post,1)-dint) (pts(ant,1)-dint) bkg(1)];
          else
            all_pts(j,:) = [(pts(post,1)-dint) (pts(ant,1)) bkg(1)];
          end
        end

        names{i,1} = files(i).name(1:find(files(i).name=='_')-1);
        group = ones(size(all_pts,1), 1)*i;
        all_vals = [all_vals; [all_pts group]];
      end

      bkg = find(ismember(names, 'n2-par2RNAi'));
      bkg = mean(all_vals(all_vals(:,end)==bkg, :));

      all_vals(:,1:2) = bsxfun(@minus, all_vals(:,1:2), bkg(1,1:2));

      ref = find(ismember(names, 'n2'));
      avg = mean(all_vals(all_vals(:,end)==ref, :));

      figure;

      subplot(2,3,1);
      boxplot(100*(all_vals(:,1)./all_vals(:,3))*(avg(3)/avg(1)), all_vals(:,end), 'labels', names);
      title('Post membrane over bkg');

      subplot(2,3,2);
      boxplot(100*(all_vals(:,2)./all_vals(:,3))*(avg(3)/avg(2)), all_vals(:,end), 'labels', names);
      title('Ant membrane over bkg');

      [H,P] = myttest(all_vals(:,1)./all_vals(:,3), all_vals(:,end))
      [H,P] = myttest(all_vals(:,2)./all_vals(:,3), all_vals(:,end))

      files = dir('signal_quantif/*.txt');
      nfiles = length(files);
      names = cell(nfiles, 1);

      objective = 'Mean';
      %objective = 'Max';

      all_vals = NaN(0,4);
      for i=1:nfiles
        tmp = importdata(['signal_quantif/' files(i).name]);
        tmp = [tmp.data(:, ismember(tmp.colheaders, 'Slice')), ...
               tmp.data(:, ismember(tmp.colheaders, objective)), ...
               tmp.data(:, ismember(tmp.colheaders, 'BX')), ...
               tmp.data(:, ismember(tmp.colheaders, 'Width'))];

        slices = unique(tmp(:,1)).';
        all_pts = NaN(length(slices), 3);

        for j=slices
          pts = tmp(tmp(:,1)==j, 2:end);
          [junk, ant] = min(pts(:,2));
          [junk, post] = max(pts(:,2)+pts(:,3));

          bkg = [1:size(pts,1)];
          bkg = bkg(bkg ~= ant & bkg ~= post);
          bkg = sort(pts(bkg, 1));

          dint = diff(bkg);

          all_pts(j,:) = [(pts(post,1)-dint) (pts(ant,1)) bkg(1)];
        end

        names{i,1} = files(i).name(1:find(files(i).name=='_')-1);
        group = ones(size(all_pts,1), 1)*i;
        all_vals = [all_vals; [all_pts group]];
      end

      %bkg = find(ismember(names, 'n2-par2RNAi'));
      bkg = mean(all_vals(:, 2));

      %all_vals(:,1:2) = bsxfun(@minus, all_vals(:,1:2), bkg(1,1:2));
      all_vals(:,[1 3]) = all_vals(:,[1 3]) - bkg;

      ref = find(ismember(names, '1054'));
      avg = mean(all_vals(all_vals(:,end)==ref, :));

      subplot(2,3,3);
      boxplot(100*(all_vals(:,1)./all_vals(:,3))*(avg(3)/avg(1)), all_vals(:,end), 'labels', names);
      title('GFP post membrane over bkg');

      [H,P] = myttest(all_vals(:,1)./all_vals(:,3), all_vals(:,end))

      %figure;subplot(2,3,4);
      %{
      all_vals = NaN(0,4);
      for i=1:nfiles
        tmp = importdata(['signal_quantif/' files(i).name]);
        %good_vals = ismember(tmp.colheaders, 'Mean') | ismember(tmp.colheaders, 'Median');
        good_vals = ismember(tmp.colheaders, objective);

        tmp = tmp.data(:, good_vals);
        tmp = [tmp(1:2:end,:) - tmp(2:2:end,:) tmp(1:2:end,:) tmp(2:2:end,:)];
        names{i,1} = files(i).name(1:find(files(i).name=='_')-1);
        group = ones(size(tmp,1), 1)*i;
        all_vals = [all_vals; [tmp group]];
      end

      ref = find(ismember(names, '1054'));
      avg = mean(all_vals(all_vals(:,end)==ref, :));

      subplot(2,3,4);
      boxplot(100*all_vals(:,1)/avg(1), all_vals(:,end), 'labels', names);
      title('Membrane minus background');
      subplot(2,3,5);
      boxplot(100*(all_vals(:,2) ./ all_vals(:,3))*(avg(3)/avg(2)), all_vals(:,end), 'labels', names);
      title('Membrane over bkg');
      subplot(2,3,6);
      boxplot(100*all_vals(:,2)/avg(2), all_vals(:,end), 'labels', names);
      title('Membrane signal');
      %}

      %scores = NaN(nfiles);
      %for i=1:nfiles
      %  for j=i+1:nfiles
      %    [junk, scores(i,j)] =  ttest2(all_vals(all_vals(:,end)==i,1), all_vals(all_vals(:,end)==j,1));
      %    [junk, scores(j,i)] =  ttest2(all_vals(all_vals(:,end)==i,2), all_vals(all_vals(:,end)==j,2));
      %  end
      %end
      %[H,P] = myttest(all_vals(:,1), all_vals(:,end))
      %[H,P] = myttest(all_vals(:,2)./all_vals(:,3), all_vals(:,end))
      %[H,P] = myttest(all_vals(:,2), all_vals(:,end))

      targets = {'good_24.txt'; 'good_20.txt'; 'good_13.txt'; 'good_c27d91.txt'; 'good_ani2.txt'};
      data = load('data_expansion.mat');
      [junk, indxs] = ismember(targets, data.all_data(:,1));
      all_data = data.all_data(indxs,2);

      all_sizes = NaN(0,4);
      for i=1:length(targets)
        %files = textread(all_files{i}, '%s');
        %tmp_size = NaN(length(files), 3);
        %for j=1:length(files)
        %  load(files{j});
        %  tmp_size(j,:) = 2*mymovie.metadata.axes_length_3d.';
        %end
        tmp_size = 2*all_data{i}(:,3:5);
        group = ones(size(tmp_size,1), 1)*i;
        all_sizes = [all_sizes; [tmp_size group]];
      end

      mins = floor(min(all_sizes));
      maxs = ceil(max(all_sizes));
      nbins = 32;

      edges = NaN(nbins, 3);
      counts = NaN(nbins, length(targets), 3);

      for i=1:3
        edges(:, i) = linspace(mins(i), maxs(i), nbins).';
        for j=1:length(targets)
          counts(:,j,i) = histc(all_sizes(all_sizes(:,end)==j, i), edges(:,i));
        end
      end

      ymax = max(squeeze(max(counts)));
      for i=1:3
        figure;
        for j=1:length(targets)
          subplot(2,3,j);
          bar(edges(:,i), counts(:,j,i), 'histc');
          ylim([0 ymax(i)]);
          legend(strrep(targets{j}, '_', '\_'));
        end
        mtit(num2str(i));
      end

      keyboard

    case 2
      load('1056-24-040311_1_');
      nimg = 30;
      opts.verbosity = 3;

      time = get_manual_timing(mymovie, opts);
      time = time(:);

      opts.recompute = true;
      dp_dic(mymovie, nimg, opts);
      dp_markers(mymovie, nimg, opts);

      opts.segment = false;
      cortical_signal(mymovie, nimg, opts);

      opts.recompute = false;
      opts.verbosity = 0;

      [domain, ruffles, theta] = gather_quantification(mymovie, opts);
      [fraction, max_width, cell_width, raw_domain, pos, path_center] = domain_expansion(mymovie, opts);

%      center = (size(raw_domain,2)-1)/2+1;

      center = find(theta == 0);
      pos_indx = [center:50:size(domain,2)];
      tmp = [center:-50:1];
      pos_indx = [tmp(end:-1:2) pos_indx];

      t_init = find(~isnan(fraction), 1, 'first');
      t_indx = unique([fliplr([t_init:-20:1]) t_init:20:size(domain,1)]);
      t_pos = ([1:size(domain,1)]-t_init)*10;

      %hold off;
      figure;
      imagesc(domain);
      hold on;
      plot(path_center, 1:length(path_center), 'k');
      plot(path_center(1:length(fraction)) + fraction*max_width*2, 1:length(fraction), 'k');
      plot(path_center(1:length(fraction)) - fraction*max_width*2, 1:length(fraction), 'k');
      plot(repmat([1 size(domain,2)], 3, 1).', time(:,[1 1]).', 'k')
      plot([1 size(domain,2)], [nimg nimg], 'w')
      title('1056-c27d91-040712\_3\_');
      set(gca, 'XTick', pos_indx, 'XTickLabel', theta(pos_indx));
      set(gca, 'YTick', t_indx, 'YTickLabel', t_pos(t_indx));
      colormap(blueredmap)

      width = (size(raw_domain, 2) - 1) / 2;
      pos_indx = [0:50:width];
      path = fraction*max_width/opts_expansion.quantification.resolution;

      t_indx = unique([fliplr([t_init:-20:1]) t_init:20:size(raw_domain,1)]);
      t_pos = ([1:size(raw_domain,1)]-t_init)*10;

      figure;imagesc(raw_domain)
      hold on;
      plot(path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      plot(-path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      plot(repmat([1 size(raw_domain,2)], 3, 1).', time(:,[1 1]).', 'k')
      plot([1 size(raw_domain,2)], [nimg nimg], 'w')
      set(gca, 'XTick', [-pos_indx([end:-1:2]) pos_indx]+width+1, 'XTickLabel', pos([-pos_indx([end:-1:2]) pos_indx]+width+1));
      set(gca, 'YTick', t_indx, 'YTickLabel', t_pos(t_indx));
      title('1056-c27d91-040712\_3\_');
      colormap(blueredmap)

      keyboard

    case 3
      t0 = 35; % corr_offset + yindx in the next find_kymograph

      advection_params = getfield(load('cyto_flow.mat', 'flow3'), 'flow3');
      figure;imagesc([-advection_params advection_params(:,[end-1:-1:1])])
      pos_indx = fliplr([size(advection_params,2):-25:1]);
      pos_indx = [pos_indx(1:end-1) pos_indx(end):25:size(advection_params,2)*2-1];
      pos = [1:2*size(advection_params,2)-1]-size(advection_params,2);

      t_indx = unique([fliplr([t0*10:-200:1]) t0*10:200:size(advection_params,1)]);

      set(gca, 'XTick', pos_indx, 'XTickLabel', pos(pos_indx))
      set(gca, 'YTick', t_indx, 'YTickLabel', t_indx - t0*10)
      colormap(redgreenmap)
      colorbar

      %figure;find_kymograph('1056-24-all.mat', 'config_fitting', 'fit_kymo', 'config_modeling', 'goehring', 'init_noise', 0, 'display', true, 'aligning_type', 'best');

      figure;find_kymograph('1056-24-all.mat', 'config_fitting', 'fit_kymo', 'config_modeling', 'custom_flow', 'init_noise', 0, 'display', true, 'aligning_type', 'best');

    case 6

      load('data_expansion.mat');
      all_times = NaN(0,5+length(thresh));

      for i=1:size(all_data, 1)
        files = all_data{i,3};
        nfiles = size(files, 1);
        tmp_time = NaN(nfiles, size(all_times,2)-1);
        for j=1:nfiles
          tmp_indx = NaN(size(thresh));
          fraction = files{j,2};
          for t=1:length(thresh)
            tmp_indx(t) = find(fraction>=thresh(t), 1, 'first');
          end
          tmp_time(j,:) = [length(fraction) tmp_indx all_data{i,2}(j,6:8)];
        end
        group = ones(size(tmp_time,1), 1)*i;
        all_times = [all_times; [tmp_time group]];
      end

      figure;
      subplot(2,3,1)
      boxplot(all_times(:,1), all_times(:,end), 'label', all_data(:,1));
      title('Recording duration')
      [H,P] = myttest(all_times(:,1), all_times(:,end))

      subplot(2,3,2)
      boxplot(all_times(:,6)-all_times(:,4), all_times(:,end), 'label', all_data(:,1));
      title('PC to CK')
      [H,P] = myttest(all_times(:,6)-all_times(:,4), all_times(:,end))

      subplot(2,3,3)
      boxplot(all_times(:,5)-all_times(:,4), all_times(:,end), 'label', all_data(:,1));
      title('PC to PNM')
      [H,P] = myttest(all_times(:,5)-all_times(:,4), all_times(:,end))
      
      subplot(2,3,4)
      boxplot(all_times(:,3)-all_times(:,2), all_times(:,end), 'label', all_data(:,1));
      title('Polarity establishment')
      [H,P] = myttest(all_times(:,3)-all_times(:,2), all_times(:,end))

      subplot(2,3,5)
      boxplot((all_times(:,3)-all_times(:,2)) ./ (all_times(:,5)-all_times(:,4)), all_times(:,end), 'label', all_data(:,1));
      title('Polarity establishment VS PC to PNM')
      [H,P] = myttest((all_times(:,3)-all_times(:,2)) ./ (all_times(:,5)-all_times(:,4)), all_times(:,end))

      subplot(2,3,6)
      %scatter((all_times(:,3)-all_times(:,2)), (all_times(:,5)-all_times(:,4)));
      myregress([ones(size(all_times,1), 1) (all_times(:,3)-all_times(:,2))], (all_times(:,5)-all_times(:,4)));
      title('Polarity establishment VS PC to PNM')

      %find_kymograph('1056-temps-all.mat', 'config_modeling', 'custom_flow', 'config_fitting', 'fit_kymo', 'display', true);

      keyboard

    case 6.1

      colors = [179 205 227; ...
                204 235 197; ...
                251 180 174; ...
                55 126 184; ...
                77 175 74; ...
                228 26 28]/255;

      all_data = cell(3, 3);

      data = load('data_expansion.mat');

      targets = {'good_13.txt';'good_20.txt';'good_24.txt'};

      for f=1:size(data.all_data, 1)
        files = data.all_data{f,3};
        nfiles = size(files,1);

        all_maint = cell(nfiles, 1);
        centers = NaN(1, nfiles);
        all_params = NaN(nfiles, 6);

        for i=1:nfiles
          all_maint{i} = files{i,1}.';
          centers(i) = round((data.all_data{f,2}(i,1)/2)/opts_expansion.quantification.resolution);

          tmp_indx = NaN(size(thresh));
          fraction = files{i,2};
          for t=1:length(thresh)
            tmp_indx(t) = find(fraction>=thresh(t), 1, 'first');
          end

          all_params(i,:) = [data.all_data{f,2}(i,1:5) (tmp_indx(2)-tmp_indx(1))*10];
          disp([num2str(i) '/' num2str(nfiles)]);
        end
        [stack, shift] = stack_images(all_maint, centers, 0.5);
        center = centers(1) + shift(1);

        all_data{f,1} = stack;
        all_data{f,2} = all_params;
        all_data{f,3} = center;
      end

      [junk, indxs] = ismember(targets, data.all_data(:,1));
      temperatures = data.temperatures(indxs)
      all_data = all_data(indxs,:);
      pos = cell(1,3);

      all_times = NaN(0,2);

      figure;hold on
      for i=1:3
        boxplot(all_data{i,2}(:,6), 'colors',colors(i, :), 'position', temperatures(i), 'width', 2);
        all_times = [all_times; [all_data{i,2}(:,6) ones(size(all_data{i,2}),1)*i]];
      end
      avgs = mymean(all_times(:,1), 1, all_times(:,2));
      plot(temperatures, avgs, '-ok');

      title('Polarity establishment')
      ylim([0 max(all_times(:,1)+10)]);
      xlim([10 26]);
      set(gca, 'XTick', [10:26], 'XTickLabel', [10:26]);
      [H,P] = myttest(all_times(:,1), all_times(:,end))

      figure;
      for i=1:3
        stack = all_data{i,1};
        center = all_data{i,3};

        [m24, s24] = mymean(stack, 3);
        max_val = mean(m24(1:ntails));
        min_val = mean(m24(end-ntails+1:end));
        norm_factor = 1 / (max_val - min_val);

        m24 = (m24 - min_val) *norm_factor;
        s24 = s24 * norm_factor;
        stack = (stack - min_val) * norm_factor;

        all_data{i,1} = stack;

        pos{i} = ([1:size(stack,1)]-center)*opts_expansion.quantification.resolution;;

        subplot(2,3,i)
        hold on;
        plot(pos{i}, squeeze(stack),  'Color', colors(i, :));
        plot(pos{i}, m24+s24, 'Color', colors(i+3,:));
        plot(pos{i}, m24-s24, 'Color', colors(i+3,:));
        plot(pos{i}, m24, 'Color', colors(i+3,:), 'LineWidth', 2);
        ylim([0 1.2]);
        xlim([-20 20]);
        xlabel('Distance to the boundary (µm)');
        ylabel('GFP intensity (a.u.)')
        switch i
          case 1
            title('13°C')
          case 2
            title('20°C')
          case 3
            title('24°C')
        end

        subplot(2,3,4);hold on
        plot(pos{i}, squeeze(stack),  'Color', colors(i, :));
        %plot(pos{i}, m24+s24,  'Color', colors(i, :));
        %plot(pos{i}, m24-s24,  'Color', colors(i, :));
        plot(pos{i}, m24, 'Color', colors(i+3, :), 'LineWidth', 2);
      end
      ylim([0 1.2]);
      xlim([-20 20]);

      [all_stacks, all_shift] = stack_images(all_data(:,1), cat(1,all_data{:,3}), 0.5);
      [m_all, s_all] = mymean(all_stacks, 3);
      pos_all = ([1:size(all_stacks,1)]-(all_data{1,3}+all_shift(1)))*opts_expansion.quantification.resolution;;

      subplot(2,3,4);hold on;
%      plot(pos_all, m_all+s_all, 'Color', [83 83 83]/255);
%      plot(pos_all, m_all-s_all, 'Color', [83 83 83]/255);
      plot(pos_all, m_all, 'Color', [83 83 83]/255, 'LineWidth', 2);
      xlabel('Distance to the boundary (µm)');
      ylabel('GFP intensity (a.u.)')
      title('Average profiles')
      ylim([0 1.2]);
      xlim([-20 20]);

      subplot(2,3,5)
      hold on;
      plot(pos_all, squeeze(all_stacks), 'Color', [189 189 189]/255);
      plot(pos_all, m_all+s_all, 'Color', [83 83 83]/255);
      plot(pos_all, m_all-s_all, 'Color', [83 83 83]/255);
      plot(pos_all, m_all, 'Color', [83 83 83]/255, 'LineWidth', 2);
      xlabel('Distance to the boundary (µm)');
      ylabel('GFP intensity (a.u.)')
      title('Overall average profile')
      ylim([0 1.2]);
      xlim([-20 20]);

      disp([adaptive_neyman(squeeze(all_data{1,1}), squeeze(all_data{2,1})), ...
            adaptive_neyman(squeeze(all_data{1,1}), squeeze(all_data{3,1})), ...
            adaptive_neyman(squeeze(all_data{2,1}), squeeze(all_data{3,1}))]);
      
      figure;
      for i=1:3
        all_params = all_data{i,2};
        intercept = ones(size(all_params,1), 1);
        [means, stds] = mymean([2*all_params(:,3) all_params(:,1)./all_params(:,2)]);

        subplot(2,3,i);hold on
        ylim([0 1]);
        xlim([30 80]);
        errorbarxy(means(1), means(2), stds(1), stds(2));hold on;
        myregress([intercept 2*all_params(:,3)], all_params(:,1)./all_params(:,2), colors(i, :));
        title(['N=' num2str(length(intercept))])
        xlabel('Egg length (µm)');
        ylabel('Fraction (a.u.)')

        subplot(2,3,4);hold on
        ylim([0 1]);
        xlim([30 80]);
        errorbarxy(means(1), means(2), stds(1), stds(2)); hold on
        myregress([intercept 2*all_params(:,3)], all_params(:,1)./all_params(:,2), colors(i, :));
      end
      title('13, 20, 24')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')
      ylim([0 1]);
      xlim([30 80]);

      all_params = cat(1, all_data{:,2});
      intercept = ones(size(all_params,1), 1);

      subplot(2,3,5);hold on
      ylim([0 1]);
      xlim([30 80]);
      myregress([intercept 2*all_params(:,3)], all_params(:,1)./all_params(:,2));
      title('All together')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')

      mtit('Membrane length fraction')

      keyboard

    case 6.2

      opts = get_struct('modeling');
      opts = load_parameters(opts, 'goehring.txt');
      opts = load_parameters(opts, 'custom_flow.txt');

      npos = opts.nparticles;

      flow = opts.advection_params;
      if (size(flow, 1) ~= npos)
        [X, Y] = meshgrid([1:size(flow, 2)], 1+([0:npos-1]*(size(flow, 1)-1)/(npos-1)).');
        flow = bilinear_mex(flow, X, Y, [2 2]);
      end

      orig = opts;

      %scalings = [0.25:0.25:1.75];
      %scalings = [0.1 0.25 0.5 1 2 4 10];
      %scalings = [0.25 0.33 0.5 0.75 1 1.5 2 3 4];
      % because: exp(-(E/kB)*((1/(13+C2K)) - (1/(20+C2K)))) = 0.5329
      % and      exp(-(E/kB)*((1/(24+C2K)) - (1/(20+C2K)))) = 1.4139
      scalings = [0.1:0.1:1.5];

      colors = redbluemap(length(scalings));
      colors(((length(scalings)-1)/2)+1,:) = [0 0 0];

      labels = cellstr(num2str(scalings.'));

      opts.max_iter = opts.max_iter * 5;


      x0 = opts.init_func(opts, false);
      ml_params = [opts.diffusion_params; ...
                      opts.reaction_params];

      opts_expansion = load_parameters(get_struct('ASSET'), 'domain_expansion.txt');

      nsimul = length(scalings);
      all_maint = cell(nsimul, 1);
      centers = NaN(1, nsimul);
      all_params = NaN(nsimul, 6);

      titles = {'D_{A/P}', 'k_{on}', 'k_{off}', 'k_{AP/PA}', '\alpha / \beta', '\nu'};
      colors = redbluemap(nsimul);

      for i=1:6
        if (mod(i-1, 3)==0)
          figure;hold on;
        end
        for j=1:length(scalings)
          flow_scale = 1;
          switch i
            case 1
              opts.diffusion_params = orig.diffusion_params * scalings(j);
            case 2
              opts.diffusion_params = orig.diffusion_params;
              opts.reaction_params(1,:) = orig.reaction_params(1,:) * scalings(j);
            case 3
              opts.reaction_params(1,:) = orig.reaction_params(1,:);
              opts.reaction_params(2,:) = orig.reaction_params(2,:) * scalings(j);
            case 4
              opts.reaction_params(2,:) = orig.reaction_params(2,:);
              opts.reaction_params(3,:) = orig.reaction_params(3,:) * scalings(j);
            case 5
              opts.reaction_params(3,:) = orig.reaction_params(3,:);
              opts.reaction_params(4,:) = orig.reaction_params(4,:) * scalings(j);
            case 6
              opts.reaction_params(4,:) = orig.reaction_params(4,:);
              flow_scale = scalings(j);
          end

          x0 = opts.init_func(opts, false);
          ml_params = [opts.diffusion_params; ...
                        opts.reaction_params];

          [simul, t_pos] = simulate_model_mix(x0, ml_params, opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow_scale*flow, opts.user_data, opts.max_iter);

          domain = simul(npos+1:end, :).';
          domain = [domain domain(:,end-1:-1:1)];

          [profile, center, max_width, cell_width, path] = get_profile(domain, nframes);
          timing = (find(path>=thresh(2), 1, 'first') - (find(path>=thresh(1), 1, 'first')))*10;

          all_maint{j} = profile.';
          centers(j) = center;
          all_params(j,:) = [2*max_width, cell_width, opts.axes_length(:).' timing];
          all_params(j,1:2) = all_params(j,1:2)*opts.x_step/opts_expansion.quantification.resolution;
          
          %subplot(4,3,mod(i-1,3)+10);
          %hold on;
          %plot(path, 'Color', colors(j,:))

          disp([num2str(j) '/' num2str(nsimul) ':' num2str(i) '/6']);
        end
        %if (i==3)
        %keyboard
        %end

        goods = (centers~=npos);
        temp_colors = colors(goods,:);

        centers = centers(goods);

        [stack, shift] = stack_images(all_maint(goods), centers, 0.5);
        center = centers(1) + shift(1);
        [m_simul, s_simul] = mymean(stack, 3);
        pos_simul = ([1:size(stack,1)]-center)*opts.x_step;

        intercept = ones(size(all_params,1), 1);

        hsub = subplot(3,3,mod(i-1,3)+1);
        set(hsub,'ColorOrder', temp_colors);
        hold on;
      

        plot(pos_simul, squeeze(stack));
        plot(pos_simul, m_simul+s_simul, 'k');
        plot(pos_simul, m_simul-s_simul, 'k');
        plot(pos_simul, m_simul, 'k', 'LineWidth', 2);
        xlabel('Distance to the boundary (µm)');
        ylabel('Intensity rescaled (a.u.)');
        title('Average profile')
        ylim([0 1]);
        xlim([-20 20]);
        title(titles(i))
        legend(labels(goods));

        subplot(3,3,mod(i-1,3)+4);hold on
        myregress([intercept(goods) all_params(goods,3)], all_params(goods,1)./all_params(goods,2), temp_colors);
        title('Membrane fraction')
        xlabel('Egg length (µm)');
        ylabel('Fraction (a.u.)')
        ylim([0 1]);

        subplot(3,3,mod(i-1,3)+7);hold on
        plot(scalings(goods), all_params(goods,end))

        %keyboard

      end

      keyboard

      return;


      opts = get_struct('modeling');
      opts = load_parameters(opts, 'goehring.txt');
      opts = load_parameters(opts, 'custom_flow.txt');
      opts = load_parameters(opts, 'maintenance.txt');

      %% Dividing both equations by a constant is equal to time rescaling, which
      %% means that nothing changes if we simulate long enough, hence less parameters 
      %% to check
      opts.diffusion_params = opts.diffusion_params ./ opts.reaction_params(3, :);
      opts.reaction_params(1:3,:) = bsxfun(@rdivide, opts.reaction_params(1:3,:), opts.reaction_params(3, :));

      flow = opts.advection_params;
      if (size(flow, 1) ~= opts.nparticles)
        [X, Y] = meshgrid([1:size(flow, 2)], 1+([0:opts.nparticles-1]*(size(flow, 1)-1)/(opts.nparticles-1)).');
        flow = bilinear_mex(flow, X, Y, [2 2]);
      end

      orig = opts;
      %scalings = [0.25:0.25:1.75];
      %scalings = [0.1 0.25 0.5 1 2 4 10];
      scalings = [0.25 0.33 0.5 0.75 1 1.5 2 3 4];

      colors = redbluemap(length(scalings));
      colors(((length(scalings)-1)/2)+1,:) = [0 0 0];

      labels = cellstr(num2str(scalings.'));

      opts.max_iter = opts.max_iter * 3;
      titles = {'D_{A/P}', 'k_{on}', 'k_{off}', '\alpha / \beta'};

      figure;hold on;
      for i=1:4
      %for i=1:1
        subplot(2,2,i);hold on;
        title(num2str(i))
        for j=1:length(scalings)
          switch i
            case 1
              opts.diffusion_params = orig.diffusion_params * scalings(j);
            case 2
              opts.diffusion_params = orig.diffusion_params;
              opts.reaction_params(1,:) = orig.reaction_params(1,:) * scalings(j);
            case 3
              opts.reaction_params(1,:) = orig.reaction_params(1,:);
              opts.reaction_params(2,:) = orig.reaction_params(2,:) * scalings(j);
            case 4
              opts.reaction_params(2,:) = orig.reaction_params(2,:);
              opts.reaction_params(4,:) = orig.reaction_params(4,:) * scalings(j);
          end

          x0 = opts.init_func(opts, false);
          ml_params = [opts.diffusion_params; ...
                        opts.reaction_params];

          [simul, t_pos] = simulate_model_mix(x0, ml_params, opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow, opts.user_data, opts.max_iter);

          domain = simul((end/2)+1:end, end).';
          npts = length(domain);
          %domain = [domain domain(end-1:-1:1)];
          pos = opts.reaction_params(end) * [npts-1:-1:0] / (npts-1);

          plot(pos, 100*domain/(max(domain)), 'Color', colors(j,:));

          disp([num2str(j) '/' num2str(i)])
        end
        xlim([0 opts.reaction_params(end)])
        ylabel('Intensity rescaled (a.u.)');
        xlabel('Distance to pole (µm)')
        title(titles(i));
        legend(labels);
      end

      keyboard

    case 7
      %vals = group_ml_results('LatestFits/adr-kymo-*_evol.dat', {'parameter_set';'fit_flow';'fit_model'}, {'type', '1056-temps-all'; 'fitting_type', 'cmaes'; 'aligning_type', 'fitting';'normalize_smooth', true; 'rescale_length_only', true});
      vals = group_ml_results('LatestFits/adr-kymo-*_evol.dat', {'parameter_set';'fit_flow';'fit_model';'scale_flow'}, {'type', '1056-temps-all'; 'aligning_type', 'fitting';'normalize_smooth', true; 'rescale_length_only', true; 'scale_each_egg', true});

      data = load('1056-temps-all.mat');
      npts = sum(cellfun(@(x)(size(x,1)), data.ground_truth));
      %load('data_fitting.mat');
      %good = ismember(all_data(:,1), 'averages');
      %values = all_data(good, :);
      %good = ismember(values{1,3}, '1056-temps-all.mat');
      %npts = values{1,2}(good, end);

      param_set = NaN(size(vals,1), 4);
      score = NaN(size(param_set, 1), 1);
      nparams = score;
      params = cell(size(score));
      rel_params = params;
      orig_p = params;

      orig_v = vals;
      vals = extract_model_parameters(vals, false);

      for i=1:size(vals,1)
        param_set(i,:) = [vals{i,1}{1}.parameter_set vals{i,1}{1}.fit_flow vals{i,1}{1}.fit_model vals{i,1}{1}.scale_flow];
        best = Inf;
        indx = 0;
        for j=1:size(vals{i,2}, 1)

          %[v, m] = deal(vals{i,2}{j,2}, vals{i,1}{1});
          %p = v(end).params;
          %ngroups = length(m.temperature);
          %p = p .* [m.rescale_factor m.offset_scaling*ones(1,ngroups) ones(1,length(p) - 4 - ngroups)];

          %p = [vals{i,2}{j,2}.params.rate vals{i,2}{j,2}.params.offset vals{i,2}{j,2}.params.energy vals{i,2}{j,2}.params.viscosity vals{i,2}{j,2}.params.flow vals{i,2}{j,2}.params.sigma]

          %disp(vals{i,2}{j,1})
          %disp(vals{i,2}{j,2}.score)

          %old_score = vals{i,2}{j,2}.score;
          %new_score = check_parameters([vals{i,1}{1}.type '.mat'], vals{i,1}{1}, 'config_modeling', 'custom_flow', 'init_pos', p, 'start_with_best', false);

          %if (abs(old_score - new_score) > 2)
          %  beep;keyboard
          %end
  
          %find_kymograph([m.type '.mat'], m, 'config_modeling', 'custom_flow', 'init_pos', p, 'init_noise', 0, 'display', true);
          %find_kymograph([m.type '.mat'], m, 'config_modeling', 'custom_flow', 'init_pos', p2, 'init_noise', 0, 'display', true);

          %keyboard

          if (best > vals{i,2}{j,2}.score)
            %disp([param_set(i,:) length(vals{i,2}{j,2}(end).params)])
          %  if (~(param_set(i,1)==24 && length(vals{i,2}{j,2}(end).params)>12))
              best = vals{i,2}{j,2}.score;
              indx = j;
          %  end
          end
        end

        if (indx ~= 0)
          score(i) = vals{i,2}{indx,2}.score;
          %p = vals{i,2}{indx,2}(end).params;
          %p(1:4) = abs(p(1:4)) .* vals{i,1}{1}.rescale_factor;
          %orig_p{i} = p;
          %p([1 3]) = p([1 3]) .* p([4 2]);
          %params{i} = p;
          %nparams(i) = numel(params{i});

          params{i} = [vals{i,2}{indx,2}.params.rate vals{i,2}{indx,2}.params.offset vals{i,2}{indx,2}.params.energy vals{i,2}{indx,2}.params.viscosity vals{i,2}{indx,2}.params.flow vals{i,2}{indx,2}.params.sigma];
          try
          rel_params{i} = [vals{i,2}{indx,2}.params.effective_value.viscosity; vals{i,2}{indx,2}.params.effective_value.rate; vals{i,2}{indx,2}.params.effective_value.flow];
          catch
            beep;keyboard
          end
          nparams(i) = vals{i,2}{indx,2}.params.nparams;
        else
          keyboard
        end
      end

      [param_set, indx] = sortrows(param_set);
      score = score(indx);
      nparams = nparams(indx);
      params = params(indx);
      rel_params = rel_params(indx);

      aic = 2*(nparams + score) + 2*nparams.*(nparams+1)./(npts-nparams-1);

      temps = vals{1,2}{1,2}.params.temperature;

      colors = [0 0 0; ...
                37 37 37; ...
                82 82 82; ...
                150 150 150; ...
                217 217 217]/255;

      for i=1:length(aic)
        if (mod(i, 6) == 1)
          figure;
        end
        subplot(2,3,mod(i-1,6)+1);hold on;
        vals = rel_params{i};
        vals = bsxfun(@rdivide, vals, vals(:,2));

        for j=2:size(vals,1)-1
          plot(temps, vals(j,:), '-s', 'Color', colors(j-1,:));
        end
        plot(temps, vals(end,:), '-v', 'Color', colors(end-1,:));
        plot(temps, vals(1,:), '-o', 'Color', colors(end,:));
        ylim([0 2.5]);
        xlim([10 26])

        title([num2str(param_set(i,:)) ':' num2str(aic(i))]);
      end

      keyboard

    case 7.1
      vals = group_ml_results('LatestFits/adr-kymo-*_evol.dat', {'init_noise'}, {'type', 'simulation'});
      colors = [37 82 115 150] / 255;

      has_drawn_init_pos = false;

      is_log = [0 1 1 1 1 0];

      noise = NaN(size(vals,1),1);
      all_best = cell(size(noise));
      max_iter = 0;
      init_pos = [];

      for i=1:length(vals)
        noise(i) = vals{i,1}{1}.init_noise;
        best_vals = NaN(0, 6);
        for j=1:size(vals{i,2}, 1)
          best_vals(end+1,:) = [vals{i,2}{j,2}(end).score vals{i,2}{j,2}(end).params];
          niter = size(vals{i,2}{j,2}(end).evolution{1},1);
          if (niter>max_iter)
            max_iter = niter;
          end
        end
        if (isempty(init_pos) && ~isempty(vals{i,1}{1}.init_pos))
          init_pos = vals{i,1}{1}.init_pos;
        end
        all_best{i} = best_vals;
      end
      [junk, indx] = sort(noise, 'descend');
      vals = vals(indx,:);
      all_best = all_best(indx);
      best_vals = cat(1, all_best{:});

      best_vals(:,2:5) = abs(best_vals(:,2:5));
      outliers = pcout(best_vals(:,1), [0.33 5], [0.25 1]);
      count = 1;

      init_pos = [NaN init_pos 0];
      count = 1;
      figure;hold on
      for i=1:length(vals)
        for j=1:size(vals{i,2}, 1)
          if (outliers(count))
            c = [179 179 179]/255;
          else
            c = [1 0 0] * colors(i);
          end

          pts = vals{i,2}{j,2}(end).evolution{1};
          pts(:,2:5) = abs(bsxfun(@rdivide, pts(:,2:5), vals{i,1}{1}.simulation_parameters ./ vals{i,1}{1}.rescale_factor));
          pts(:,6) = pts(:,6)*vals{i,1}{1}.offset_scaling;

          best_vals(count,2:5) = abs(bsxfun(@rdivide, best_vals(count,2:5), vals{i,1}{1}.simulation_parameters ./ vals{i,1}{1}.rescale_factor));
          best_vals(count,6) = best_vals(count,6)*vals{i,1}{1}.offset_scaling;
          init_pos(2:5) = init_pos(2:5) ./ vals{i,1}{1}.simulation_parameters;

          pos = [1:size(pts,1)] - size(pts,1);
          for n=1:size(pts,2)
            subplot(2,3,n);hold on
            plot(pos,pts(:,n), 'Color', c);
            if (n>1&&n<=5)
              set(gca,'YScale', 'log', 'YTickLabel',roundn(get(gca,'YTick'), -2));
              if (mod(n,2)==1)
                ylim([1e-1 1e1]);
              else
                ylim([1e-3 1e1]);
              end
            end
            if (i*j==1)
              plot([1 max_iter]-max_iter, [1 1]*init_pos(n), 'b');
            end
            xlim([1 max_iter]-max_iter);
          end
          count = count+1;
        end
      end

      alpha = 0.05;
      bounds = prctile(best_vals(~outliers,:), [alpha 0.5 1-alpha]*100);

      figure;
      for i=1:size(best_vals,2)

        if (is_log(i))
          bound = max(log10(max(abs(bounds(:,logical(is_log))))*1.75));
          pos = logspace(-bound, bound, 16);
        else
          med = bounds(2,i);
          bound = max(abs(bounds([1 3],i) - med)*1.75);
          pos = linspace(med - bound, med + bound, 16);
        end

        %log_pos = logspace(-0.3, 0.3, 16);
        %tmp_pos = [-Inf log_pos(2:end-1) Inf];
        tmp_pos = [-Inf pos(2:end-1) Inf];

        nhist = histc(best_vals(~outliers,i), tmp_pos);

        subplot(2,3,i);hold on
        %bar(log_pos, nhist(:,i),'histc');
        %bar(log_pos, nhist(:,i));
        bar(pos, nhist, 'histc');
        %if (i>1&&i<=5)
        if (is_log(i))
          set(gca,'XScale', 'log', 'XLim', pos([1 end]), 'XTick', diff(pos)/2+pos(1:end-1));
          set(gca, 'XTickLabel',roundn(get(gca,'XTick'), -2));
        else
          set(gca, 'XLim', pos([1 end]), 'XTick', diff(pos)/2+pos(1:end-1));
        end

        %hist(best_vals(~outliers,i), 16);
      end


      keyboard

    case 4
      rescale_size = [400 700];

      load('1056-ani2-250612_3_.mat');
      nimg = 61;

      img = imnorm(double(load_data(mymovie.dic,nimg)));
      img = mask_neighbors(img, mymovie.dic.centers(:,nimg), mymovie.dic.axes_length(:,nimg), mymovie.dic.orientations(1,nimg), mymovie.dic.neighbors(nimg), opts);
      [f, frac_width, full_width, domain, pos] = domain_expansion(mymovie, opts);
      width = (size(domain, 2) - 1) / 2;
      pos_indx = [0:50:width];
      path = f*frac_width/opts_expansion.quantification.resolution;

      yindx = find(path>0, 1, 'first');
      y_tick = unique([fliplr([yindx:-20:1]) yindx:20:size(domain,1)]);
      y_labels = (y_tick - yindx)*10;

      figure;
      imshow(realign(img,rescale_size,mymovie.dic.centers(:,nimg),mymovie.dic.orientations(1,nimg)));
      title('1056-ani2-250612\_3\_');

      figure;imagesc(domain)
      hold on;
      plot(path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      plot(-path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      set(gca, 'XTick', [-pos_indx([end:-1:2]) pos_indx]+width+1, 'XTickLabel', pos([-pos_indx([end:-1:2]) pos_indx]+width+1));
      set(gca, 'YTick', y_tick, 'YTickLabel', y_labels);
      title('1056-ani2-250612\_3\_');
      colormap(blueredmap)

      domain = permute([domain(:, 1:width+1, :), domain(:,end:-1:width+1, :)], [2 1 3]);
      noise = estimate_noise(domain);
      domain = min_max_domain(domain, path, 3*noise(:,2));
      domain = permute([domain(1:width+1, :, :); domain(end-1:-1:width+2, :, :)], [2 1 3]); 

      figure;imagesc(domain)
      hold on;
      plot(path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      plot(-path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      set(gca, 'XTick', [-pos_indx([end:-1:2]) pos_indx]+width+1, 'XTickLabel', pos([-pos_indx([end:-1:2]) pos_indx]+width+1));
      set(gca, 'YTick', y_tick, 'YTickLabel', y_labels);
      title('1056-ani2-250612\_3\_ normalized');
      colormap(blueredmap)

      load('1056-24-250511_0_.mat');
      nimg = 63;

      img = imnorm(double(load_data(mymovie.dic,nimg)));
      img = mask_neighbors(img, mymovie.dic.centers(:,nimg), mymovie.dic.axes_length(:,nimg), mymovie.dic.orientations(1,nimg), mymovie.dic.neighbors(nimg), opts);
      [f, frac_width, full_width, domain, pos] = domain_expansion(mymovie, opts);
      width = (size(domain, 2) - 1) / 2;
      pos_indx = [0:50:width];
      path = f*frac_width/opts_expansion.quantification.resolution;

      yindx = find(path>0, 1, 'first');
      y_tick = unique([fliplr([yindx:-20:1]) yindx:20:size(domain,1)]);
      y_labels = (y_tick - yindx)*10;

      figure;
      imshow(realign(img,rescale_size,mymovie.dic.centers(:,nimg),mymovie.dic.orientations(1,nimg)));
      title('1056-24-250511\_0\_');

      figure;imagesc(domain)
      hold on;
      plot(path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      plot(-path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      set(gca, 'XTick', [-pos_indx([end:-1:2]) pos_indx]+width+1, 'XTickLabel', pos([-pos_indx([end:-1:2]) pos_indx]+width+1));
      set(gca, 'YTick', y_tick, 'YTickLabel', y_labels);
      title('1056-24-250511\_0\_');
      colormap(blueredmap)

      domain = permute([domain(:, 1:width+1, :), domain(:,end:-1:width+1, :)], [2 1 3]);
      noise = estimate_noise(domain);
      domain = min_max_domain(domain, path, 3*noise(:,2));
      domain = permute([domain(1:width+1, :, :); domain(end-1:-1:width+2, :, :)], [2 1 3]); 

      figure;imagesc(domain)
      hold on;
      plot(path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      plot(-path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      set(gca, 'XTick', [-pos_indx([end:-1:2]) pos_indx]+width+1, 'XTickLabel', pos([-pos_indx([end:-1:2]) pos_indx]+width+1));
      set(gca, 'YTick', y_tick, 'YTickLabel', y_labels);
      title('1056-24-250511\_0\_ normalized');
      colormap(blueredmap)

      load('1056-c27d91-040712_3_.mat');
      nimg = 67;

      img = imnorm(double(load_data(mymovie.dic,nimg)));
      img = mask_neighbors(img, mymovie.dic.centers(:,nimg), mymovie.dic.axes_length(:,nimg), mymovie.dic.orientations(1,nimg), mymovie.dic.neighbors(nimg), opts);
      [f, frac_width, full_width, domain, pos] = domain_expansion(mymovie, opts);
      width = (size(domain, 2) - 1) / 2;
      pos_indx = [0:50:width];
      path = f*frac_width/opts_expansion.quantification.resolution;

      yindx = find(path>0, 1, 'first');
      y_tick = unique([fliplr([yindx:-20:1]) yindx:20:size(domain,1)]);
      y_labels = (y_tick - yindx)*10;

      figure;
      imshow(realign(img,rescale_size,mymovie.dic.centers(:,nimg),mymovie.dic.orientations(1,nimg)));
      title('1056-c27d91-040712\_3\_');

      figure;imagesc(domain)
      hold on;
      plot(path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      plot(-path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      set(gca, 'XTick', [-pos_indx([end:-1:2]) pos_indx]+width+1, 'XTickLabel', pos([-pos_indx([end:-1:2]) pos_indx]+width+1));
      set(gca, 'YTick', y_tick, 'YTickLabel', y_labels);
      title('1056-c27d91-040712\_3\_');
      colormap(blueredmap)

      domain = permute([domain(:, 1:width+1, :), domain(:,end:-1:width+1, :)], [2 1 3]);
      noise = estimate_noise(domain);
      domain = min_max_domain(domain, path, 3*noise(:,2));
      domain = permute([domain(1:width+1, :, :); domain(end-1:-1:width+2, :, :)], [2 1 3]); 

      figure;imagesc(domain)
      hold on;
      plot(path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      plot(-path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      set(gca, 'XTick', [-pos_indx([end:-1:2]) pos_indx]+width+1, 'XTickLabel', pos([-pos_indx([end:-1:2]) pos_indx]+width+1));
      set(gca, 'YTick', y_tick, 'YTickLabel', y_labels);
      title('1056-c27d91-040712\_3\_ normalized');
      colormap(blueredmap)

    case 4.1

      colors = [254 217 166; ...
                251 180 174; ...
                222 203 228; ...
                255 127 0; ...
                228 26 28; ...
                152 78 163]/255;

      all_data = cell(3, 3);

      data = load('data_expansion.mat');

      targets = {'good_ani2.txt';'good_24.txt';'good_c27d91.txt'};

      for f=1:size(data.all_data, 1)
        files = data.all_data{f,3};
        nfiles = size(files,1);

        all_maint = cell(nfiles, 1);
        centers = NaN(1, nfiles);
        all_params = NaN(nfiles, 5);

        for i=1:nfiles
          all_maint{i} = files{i,1}.';
          centers(i) = round((data.all_data{f,2}(i,1)/2)/opts_expansion.quantification.resolution);
          all_params(i,:) = data.all_data{f,2}(i,1:5);
          disp([num2str(i) '/' num2str(nfiles)]);
        end
        [stack, shift] = stack_images(all_maint, centers, 0.5);
        center = centers(1) + shift(1);

        all_data{f,1} = stack;
        all_data{f,2} = all_params;
        all_data{f,3} = center;
      end

      [junk, indxs] = ismember(targets, data.all_data(:,1));
      all_data = all_data(indxs,:);
      pos = cell(1,3);

      figure;
      for i=1:3
        stack = all_data{i,1};
        center = all_data{i,3};

        [m24, s24] = mymean(stack, 3);

        max_val = mean(m24(1:ntails));
        min_val = mean(m24(end-ntails+1:end));
        norm_factor = 1 / (max_val - min_val);

        m24 = (m24 - min_val) *norm_factor;
        s24 = s24 * norm_factor;
        stack = (stack - min_val) * norm_factor;

        all_data{i,1} = stack;

        pos{i} = ([1:size(stack,1)]-center)*opts_expansion.quantification.resolution;;

        subplot(2,3,i)
        hold on;
        plot(pos{i}, squeeze(stack),  'Color', colors(i, :));
        plot(pos{i}, m24+s24, 'Color', colors(i+3,:));
        plot(pos{i}, m24-s24, 'Color', colors(i+3,:));
        plot(pos{i}, m24, 'Color', colors(i+3,:), 'LineWidth', 2);
        ylim([0 1.2]);
        xlim([-20 20]);
        xlabel('Distance to the boundary (µm)');
        ylabel('GFP intensity (a.u.)')
        switch i
          case 1
            title('ani-2(RNAi)')
          case 2
            title('WT')
          case 3
            title('C27D9.1(RNAi)')
        end

        subplot(2,3,4);hold on
        plot(pos{i}, squeeze(stack),  'Color', colors(i, :));
        plot(pos{i}, m24, 'Color', colors(i+3, :), 'LineWidth', 2);
      end
      ylim([0 1.2]);
      xlim([-20 20]);

      [all_stacks, all_shift] = stack_images(all_data(:,1), cat(1,all_data{:,3}), 0.5);
      [m_all, s_all] = mymean(all_stacks, 3);
      pos_all = ([1:size(all_stacks,1)]-(all_data{1,3}+all_shift(1)))*opts_expansion.quantification.resolution;;

      subplot(2,3,4);hold on;
      %plot(pos_all, m_all, 'Color', [83 83 83]/255, 'LineWidth', 2);
      xlabel('Distance to the boundary (µm)');
      ylabel('GFP intensity (a.u.)')
      title('Average profiles')
      ylim([0 1.2]);
      xlim([-20 20]);

      subplot(2,3,5)
      hold on;
      plot(pos_all, squeeze(all_stacks), 'Color', [189 189 189]/255);
      plot(pos_all, m_all+s_all, 'Color', [83 83 83]/255);
      plot(pos_all, m_all-s_all, 'Color', [83 83 83]/255);
      plot(pos_all, m_all, 'Color', [83 83 83]/255, 'LineWidth', 2);
      xlabel('Distance to the boundary (µm)');
      ylabel('GFP intensity (a.u.)')
      title('Overall average profile')
      ylim([0 1.2]);
      xlim([-20 20]);

      disp([adaptive_neyman(squeeze(all_data{1,1}), squeeze(all_data{2,1})), ...
            adaptive_neyman(squeeze(all_data{1,1}), squeeze(all_data{3,1})), ...
            adaptive_neyman(squeeze(all_data{2,1}), squeeze(all_data{3,1}))]);

      avgs = NaN(3, 4);
      figure;
      for i=1:3
        all_params = all_data{i,2};
        intercept = ones(size(all_params,1), 1);

        [means, stds] = mymean([all_params(:,3)*2 all_params(:,1)./all_params(:,2)]);

        subplot(2,3,i);hold on
        title(['N=' num2str(length(intercept))])
        xlabel('Egg length (µm)');
        ylabel('Fraction (a.u.)')
        ylim([0 1]);
        xlim([30 80]);
        errorbarxy(means(1), means(2), stds(1), stds(2)); hold on
        myregress([intercept all_params(:,3)*2], all_params(:,1)./all_params(:,2), colors(i, :));

        subplot(2,3,4);hold on
        ylim([0 1]);
        xlim([30 80]);
        errorbarxy(means(1), means(2), stds(1), stds(2)); hold on
        myregress([intercept all_params(:,3)*2], all_params(:,1)./all_params(:,2), colors(i, :));
      end
      title('ani-2, wt, C27D9.1')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')
      ylim([0 1]);
      xlim([30 80]);

      all_params = cat(1, all_data{:,2});
      intercept = ones(size(all_params,1), 1);
      [means, stds] = mymean([all_params(:,3)*2 all_params(:,1)./all_params(:,2)]);

      subplot(2,3,5);hold on
      ylim([0 1]);
      xlim([30 80]);
      errorbarxy(means(1), means(2), stds(1), stds(2));hold on
      myregress([intercept all_params(:,3)*2], all_params(:,1)./all_params(:,2));
      title('All together')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')

      mtit('Membrane length fraction')

      [ratios, surfaces, volumes] = surface2volume(all_params(:,3:5).');
      heights = find_cap_size(all_params(:,3:5).', all_params(:,1).');
      [caps_surface, caps_volume] = spheroidal_cap(all_params(:,3:5).', heights);

      figure;

      subplot(2,3,1);
      ylim([0 1]);
      xlim([30 80]);
      myregress([intercept all_params(:,3)*2], (all_params(:,3)-heights.')./(2*all_params(:,3)), black);
      title('Major axis fraction')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')

      subplot(2,3,2);
      ylim([0 1]);
      xlim([30 80]);
      myregress([intercept all_params(:,3)*2], caps_surface./surfaces, black);
      title('Membrane surface fraction')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')

      subplot(2,3,3);
      ylim([0 1]);
      xlim([30 80]);
      myregress([intercept all_params(:,3)*2], caps_volume./volumes, black);
      title('Cell volume fraction')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')

      mtit('Alternative fractions')

      subplot(2,3,4);
      myregress([intercept all_params(:,3)*2], all_params(:,1), black);
      title('Domain size vs egg length')
      xlabel('Egg length (\mu m)');
      ylabel('Domain length (\mu m)')

      subplot(2,3,5);
      myregress([intercept all_params(:,2)], all_params(:,1), black);
      title('Domain size vs membrane length')
      xlabel('Membrane length (µm)');
      ylabel('Domain length (\mu m)')

      subplot(2,3,6);
      myregress(all_params(:,2), all_params(:,1), black);
      title('Domain size vs membrane length, no intercept')
      xlabel('Membrane length (µm)');
      ylabel('Domain length (\mu m)')

      keyboard

    case 4.2

      colors = [254 217 166; ...
                251 180 174; ...
                222 203 228; ...
                255 127 0; ...
                228 26 28; ...
                152 78 163]/255;

      all_data = cell(3, 3);

      %data = load('data_expansion.mat');
      data = load('simul_expansion.mat');
      %data.all_data = data.all_simul;

      targets = {'good_ani2.txt';'good_24.txt';'good_c27d91.txt'};

      for f=1:size(data.all_simul, 1)
        files = data.all_simul{f,3};
        nfiles = size(files,1);

        all_maint = cell(nfiles, 1);
        centers = NaN(1, nfiles);
        all_params = NaN(nfiles, 6);

        for i=1:nfiles
          all_maint{i} = files{i,1}.';
          centers(i) = round((data.all_simul{f,2}(i,1)/2)/opts_expansion.quantification.resolution);
          tmp_params = data.all_simul{f,2}(i,1:5);
          tmp_params = [tmp_params 2*tmp_params(3)];
          tmp_params(3) = ellipse_circum(tmp_params(3:5), tmp_params(2), true);

          all_params(i,:) = tmp_params;
          disp([num2str(i) '/' num2str(nfiles)]);
        end
        [stack, shift] = stack_images(all_maint, centers, 0.5);
        center = centers(1) + shift(1);

        all_data{f,1} = stack;
        all_data{f,2} = all_params;
        all_data{f,3} = center;
      end

      [junk, indxs] = ismember(targets, data.all_simul(:,1));
      all_data = all_data(indxs,:);
      pos = cell(1,3);

      figure;
      for i=1:3
        stack = all_data{i,1};
        center = all_data{i,3};

        [m24, s24] = mymean(stack, 3);

        max_val = mean(m24(1:ntails));
        min_val = mean(m24(end-ntails+1:end));
        norm_factor = 1 / (max_val - min_val);

        m24 = (m24 - min_val) *norm_factor;
        s24 = s24 * norm_factor;
        stack = (stack - min_val) * norm_factor;

        all_data{i,1} = stack;

        pos{i} = ([1:size(stack,1)]-center)*opts_expansion.quantification.resolution;;

        subplot(2,3,i)
        hold on;
        plot(pos{i}, squeeze(stack),  'Color', colors(i, :));
        plot(pos{i}, m24+s24, 'Color', colors(i+3,:));
        plot(pos{i}, m24-s24, 'Color', colors(i+3,:));
        plot(pos{i}, m24, 'Color', colors(i+3,:), 'LineWidth', 2);
        ylim([0 1.2]);
        xlim([-20 20]);
        xlabel('Distance to the boundary (µm)');
        ylabel('GFP intensity (a.u.)')
        switch i
          case 1
            title('ani-2(RNAi)')
          case 2
            title('WT')
          case 3
            title('C27D9.1(RNAi)')
        end

        subplot(2,3,4);hold on
        plot(pos{i}, squeeze(stack),  'Color', colors(i, :));
        plot(pos{i}, m24, 'Color', colors(i+3, :), 'LineWidth', 2);
      end
      ylim([0 1.2]);
      xlim([-20 20]);

      [all_stacks, all_shift] = stack_images(all_data(:,1), cat(1,all_data{:,3}), 0.5);
      [m_all, s_all] = mymean(all_stacks, 3);
      pos_all = ([1:size(all_stacks,1)]-(all_data{1,3}+all_shift(1)))*opts_expansion.quantification.resolution;;

      subplot(2,3,4);hold on;
      %plot(pos_all, m_all, 'Color', [83 83 83]/255, 'LineWidth', 2);
      xlabel('Distance to the boundary (µm)');
      ylabel('GFP intensity (a.u.)')
      title('Average profiles')
      ylim([0 1.2]);
      xlim([-20 20]);

      subplot(2,3,5)
      hold on;
      plot(pos_all, squeeze(all_stacks), 'Color', [189 189 189]/255);
      plot(pos_all, m_all+s_all, 'Color', [83 83 83]/255);
      plot(pos_all, m_all-s_all, 'Color', [83 83 83]/255);
      plot(pos_all, m_all, 'Color', [83 83 83]/255, 'LineWidth', 2);
      xlabel('Distance to the boundary (µm)');
      ylabel('GFP intensity (a.u.)')
      title('Overall average profile')
      ylim([0 1.2]);
      xlim([-20 20]);

      disp([adaptive_neyman(squeeze(all_data{1,1}), squeeze(all_data{2,1})), ...
            adaptive_neyman(squeeze(all_data{1,1}), squeeze(all_data{3,1})), ...
            adaptive_neyman(squeeze(all_data{2,1}), squeeze(all_data{3,1}))]);

      avgs = NaN(3, 4);
      figure;
      for i=1:3
        all_params = all_data{i,2};
        intercept = ones(size(all_params,1), 1);

        [means, stds] = mymean([all_params(:,6) all_params(:,1)./all_params(:,2)]);

        subplot(2,3,i);hold on
        title(['N=' num2str(length(intercept))])
        xlabel('Egg length (µm)');
        ylabel('Fraction (a.u.)')
        ylim([0 1]);
        xlim([30 80]);
        errorbarxy(means(1), means(2), stds(1), stds(2)); hold on
        myregress([intercept all_params(:,6)], all_params(:,1)./all_params(:,2), colors(i, :), 'dx');

        subplot(2,3,4);hold on
        ylim([0 1]);
        xlim([30 80]);
        errorbarxy(means(1), means(2), stds(1), stds(2)); hold on
        myregress([intercept all_params(:,6)], all_params(:,1)./all_params(:,2), colors(i, :), 'dx');
      end
      title('ani-2, wt, C27D9.1')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')
      ylim([0 1]);
      xlim([30 80]);

      all_params = cat(1, all_data{:,2});
      intercept = ones(size(all_params,1), 1);
      [means, stds] = mymean([all_params(:,6) all_params(:,1)./all_params(:,2)]);

      subplot(2,3,5);hold on
      ylim([0 1]);
      xlim([30 80]);
      errorbarxy(means(1), means(2), stds(1), stds(2));hold on
      myregress([intercept all_params(:,6)], all_params(:,1)./all_params(:,2), black, 'dx');
      title('All together')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')

      mtit('Simulation of membrane length fraction')

      [ratios, surfaces, volumes] = surface2volume(all_params(:,3:5).');
      heights = find_cap_size(all_params(:,3:5).', all_params(:,1).');
      [caps_surface, caps_volume] = spheroidal_cap(all_params(:,3:5).', heights);

      figure;

      subplot(2,3,1);
      ylim([0 1]);
      xlim([30 80]);
      myregress([intercept all_params(:,6)], (all_params(:,3)-heights.')./(2*all_params(:,3)), black, 'dx');
      title('Major axis fraction')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')

      subplot(2,3,2);
      ylim([0 1]);
      xlim([30 80]);
      myregress([intercept all_params(:,6)], caps_surface./surfaces, black, 'dx');
      title('Membrane surface fraction')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')

      subplot(2,3,3);
      ylim([0 1]);
      xlim([30 80]);
      myregress([intercept all_params(:,6)], caps_volume./volumes, black, 'dx');
      title('Cell volume fraction')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')

      mtit('Alternative fractions')

      subplot(2,3,4);
      myregress([intercept all_params(:,6)], all_params(:,1), black, 'dx');
      title('Domain size vs egg length')
      xlabel('Egg length (\mu m)');
      ylabel('Domain length (\mu m)')

      subplot(2,3,5);
      myregress([intercept all_params(:,2)], all_params(:,1), black, 'dx');
      title('Domain size vs membrane length')
      xlabel('Membrane length (µm)');
      ylabel('Domain length (\mu m)')

      subplot(2,3,6);
      myregress(all_params(:,2), all_params(:,1), black, 'dx');
      title('Domain size vs membrane length, no intercept')
      xlabel('Membrane length (µm)');
      ylabel('Domain length (\mu m)')
      ylim([30 90])
      xlim([80 180])


      mtit('Simulations')

      keyboard

    case 4.3

      colors = [179 205 227; ...
                204 235 197; ...
                251 180 174; ...
                55 126 184; ...
                77 175 74; ...
                228 26 28]/255;

      all_data = cell(3, 3);

      data = load('data_expansion.mat');
      %data = load('simul_expansion.mat');
      %data.all_data = data.all_simul;

      targets = {'good_13.txt';'good_20.txt';'good_24.txt'};

      for f=1:size(data.all_data, 1)
        files = data.all_data{f,3};
        nfiles = size(files,1);

        all_maint = cell(nfiles, 1);
        centers = NaN(1, nfiles);
        all_params = NaN(nfiles, 5);

        for i=1:nfiles
          all_maint{i} = files{i,1}.';
          centers(i) = round((data.all_data{f,2}(i,1)/2)/opts_expansion.quantification.resolution);
          all_params(i,:) = data.all_data{f,2}(i,1:5);
          disp([num2str(i) '/' num2str(nfiles)]);
        end
        [stack, shift] = stack_images(all_maint, centers, 0.5);
        center = centers(1) + shift(1);

        all_data{f,1} = stack;
        all_data{f,2} = all_params;
        all_data{f,3} = center;
      end

      [junk, indxs] = ismember(targets, data.all_data(:,1));
      all_data = all_data(indxs,:);
      pos = cell(1,3);

      figure;
      for i=1:3
        stack = all_data{i,1};
        center = all_data{i,3};

        [m24, s24] = mymean(stack, 3);

        max_val = mean(m24(1:ntails));
        min_val = mean(m24(end-ntails+1:end));
        norm_factor = 1 / (max_val - min_val);

        m24 = (m24 - min_val) *norm_factor;
        s24 = s24 * norm_factor;
        stack = (stack - min_val) * norm_factor;

        all_data{i,1} = stack;

        pos{i} = ([1:size(stack,1)]-center)*opts_expansion.quantification.resolution;;

        subplot(2,3,i)
        hold on;
        plot(pos{i}, squeeze(stack),  'Color', colors(i, :));
        plot(pos{i}, m24+s24, 'Color', colors(i+3,:));
        plot(pos{i}, m24-s24, 'Color', colors(i+3,:));
        plot(pos{i}, m24, 'Color', colors(i+3,:), 'LineWidth', 2);
        ylim([0 1.2]);
        xlim([-20 20]);
        xlabel('Distance to the boundary (µm)');
        ylabel('GFP intensity (a.u.)')
        switch i
          case 1
            title('13 °C')
          case 2
            title('20 °C')
          case 3
            title('24 °C')
        end

        subplot(2,3,4);hold on
        plot(pos{i}, squeeze(stack),  'Color', colors(i, :));
        plot(pos{i}, m24, 'Color', colors(i+3, :), 'LineWidth', 2);
      end
      ylim([0 1.2]);
      xlim([-20 20]);

      [all_stacks, all_shift] = stack_images(all_data(:,1), cat(1,all_data{:,3}), 0.5);
      [m_all, s_all] = mymean(all_stacks, 3);
      pos_all = ([1:size(all_stacks,1)]-(all_data{1,3}+all_shift(1)))*opts_expansion.quantification.resolution;;

      subplot(2,3,4);hold on;
      %plot(pos_all, m_all, 'Color', [83 83 83]/255, 'LineWidth', 2);
      xlabel('Distance to the boundary (µm)');
      ylabel('GFP intensity (a.u.)')
      title('Average profiles')
      ylim([0 1.2]);
      xlim([-20 20]);

      subplot(2,3,5)
      hold on;
      plot(pos_all, squeeze(all_stacks), 'Color', [189 189 189]/255);
      plot(pos_all, m_all+s_all, 'Color', [83 83 83]/255);
      plot(pos_all, m_all-s_all, 'Color', [83 83 83]/255);
      plot(pos_all, m_all, 'Color', [83 83 83]/255, 'LineWidth', 2);
      xlabel('Distance to the boundary (µm)');
      ylabel('GFP intensity (a.u.)')
      title('Overall average profile')
      ylim([0 1.2]);
      xlim([-20 20]);

      disp([adaptive_neyman(squeeze(all_data{1,1}), squeeze(all_data{2,1})), ...
            adaptive_neyman(squeeze(all_data{1,1}), squeeze(all_data{3,1})), ...
            adaptive_neyman(squeeze(all_data{2,1}), squeeze(all_data{3,1}))]);

      avgs = NaN(3, 4);
      figure;
      for i=1:3
        all_params = all_data{i,2};
        intercept = ones(size(all_params,1), 1);

        [means, stds] = mymean([all_params(:,3)*2 all_params(:,1)./all_params(:,2)]);

        subplot(2,3,i);hold on
        title(['N=' num2str(length(intercept))])
        xlabel('Egg length (µm)');
        ylabel('Fraction (a.u.)')
        ylim([0 1]);
        xlim([30 80]);
        errorbarxy(means(1), means(2), stds(1), stds(2)); hold on
        myregress([intercept all_params(:,3)*2], all_params(:,1)./all_params(:,2), colors(i, :));

        subplot(2,3,4);hold on
        ylim([0 1]);
        xlim([30 80]);
        errorbarxy(means(1), means(2), stds(1), stds(2)); hold on
        myregress([intercept all_params(:,3)*2], all_params(:,1)./all_params(:,2), colors(i, :));
      end
      title('13, 20, 24 °C')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')
      ylim([0 1]);
      xlim([30 80]);

      all_params = cat(1, all_data{:,2});
      intercept = ones(size(all_params,1), 1);
      [means, stds] = mymean([all_params(:,3)*2 all_params(:,1)./all_params(:,2)]);

      subplot(2,3,5);hold on
      ylim([0 1]);
      xlim([30 80]);
      errorbarxy(means(1), means(2), stds(1), stds(2));hold on
      myregress([intercept all_params(:,3)*2], all_params(:,1)./all_params(:,2));
      title('All together')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')

      mtit('Membrane length fraction')

      [ratios, surfaces, volumes] = surface2volume(all_params(:,3:5).');
      heights = find_cap_size(all_params(:,3:5).', all_params(:,1).');
      [caps_surface, caps_volume] = spheroidal_cap(all_params(:,3:5).', heights);

      figure;

      subplot(2,3,1);
      ylim([0 1]);
      xlim([30 80]);
      myregress([intercept all_params(:,3)*2], (all_params(:,3)-heights.')./(2*all_params(:,3)), black);
      title('Major axis fraction')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')

      subplot(2,3,2);
      ylim([0 1]);
      xlim([30 80]);
      myregress([intercept all_params(:,3)*2], caps_surface./surfaces, black);
      title('Membrane surface fraction')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')

      subplot(2,3,3);
      ylim([0 1]);
      xlim([30 80]);
      myregress([intercept all_params(:,3)*2], caps_volume./volumes, black);
      title('Cell volume fraction')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')

      mtit('Alternative fractions')

      subplot(2,3,4);
      myregress([intercept all_params(:,3)*2], all_params(:,1), black);
      title('Domain size vs egg length')
      xlabel('Egg length (\mu m)');
      ylabel('Domain length (\mu m)')

      subplot(2,3,5);
      myregress([intercept all_params(:,2)], all_params(:,1), black);
      title('Domain size vs membrane length')
      xlabel('Membrane length (µm)');
      ylabel('Domain length (\mu m)')

      subplot(2,3,6);
      myregress(all_params(:,2), all_params(:,1), black);
      title('Domain size vs membrane length, no intercept')
      xlabel('Membrane length (µm)');
      ylabel('Domain length (\mu m)')

      colors = [179 205 227; ...
                204 235 197; ...
                251 180 174; ...
                254 217 166; ...
                222 203 228; ...
                55 126 184; ...
                77 175 74; ...
                228 26 28; ...
                255 127 0; ...
                152 78 163]/255;



      all_data = cell(5, 3);

      data = load('data_expansion.mat');
      %data = load('simul_expansion.mat');
      %data.all_data = data.all_simul;

      targets = {'good_13.txt';'good_20.txt';'good_24.txt';'good_ani2.txt'; 'good_c27d91.txt'};

      for f=1:size(data.all_data, 1)
        files = data.all_data{f,3};
        nfiles = size(files,1);

        all_maint = cell(nfiles, 1);
        centers = NaN(1, nfiles);
        all_params = NaN(nfiles, 5);

        for i=1:nfiles
          all_maint{i} = files{i,1}.';
          centers(i) = round((data.all_data{f,2}(i,1)/2)/opts_expansion.quantification.resolution);
          all_params(i,:) = data.all_data{f,2}(i,1:5);
          disp([num2str(i) '/' num2str(nfiles)]);
        end
        [stack, shift] = stack_images(all_maint, centers, 0.5);
        center = centers(1) + shift(1);

        all_data{f,1} = stack;
        all_data{f,2} = all_params;
        all_data{f,3} = center;
      end

      [junk, indxs] = ismember(targets, data.all_data(:,1));
      all_data = all_data(indxs,:);
      pos = cell(1,5);

      avgs = NaN(3, 4);
      figure;
      for i=1:5
        all_params = all_data{i,2};
        intercept = ones(size(all_params,1), 1);

        [means, stds] = mymean([all_params(:,3)*2 all_params(:,1)./all_params(:,2)]);

        subplot(3,3,i);hold on
        title(['N=' num2str(length(intercept))])
        xlabel('Egg length (µm)');
        ylabel('Fraction (a.u.)')
        ylim([0 1]);
        xlim([30 80]);
        errorbarxy(means(1), means(2), stds(1), stds(2)); hold on
        myregress([intercept all_params(:,3)*2], all_params(:,1)./all_params(:,2), colors(i, :));

        subplot(3,3,6);hold on
        ylim([0 1]);
        xlim([30 80]);
        errorbarxy(means(1), means(2), stds(1), stds(2)); hold on
        myregress([intercept all_params(:,3)*2], all_params(:,1)./all_params(:,2), colors(i, :));
      end
      title('13, 20, 24 °C, ani2, c27d9.1')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')
      ylim([0 1]);
      xlim([30 80]);

      all_params = cat(1, all_data{:,2});
      intercept = ones(size(all_params,1), 1);
      [means, stds] = mymean([all_params(:,3)*2 all_params(:,1)./all_params(:,2)]);

      subplot(3,3,7);hold on
      ylim([0 1]);
      xlim([30 80]);
      errorbarxy(means(1), means(2), stds(1), stds(2));hold on
      myregress([intercept all_params(:,3)*2], all_params(:,1)./all_params(:,2));
      title('All together')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')

      mtit('Membrane length fraction')

    case 5

      colors = [179 205 227; ...
                204 235 197; ...
                251 180 174; ...
                254 217 166; ...
                222 203 228; ...
                179 179 179; ...
                55 126 184; ...
                77 175 74; ...
                228 26 28; ...
                255 127 0; ...
                152 78 163;
                38 38 38]/255;


      targets = {'good_13.txt';'good_20.txt';'good_24.txt';'good_ani2.txt'; 'good_c27d91.txt'; 'averages'};
      %targets = {'good_24.txt';'good_ani2.txt'; 'good_c27d91.txt'; 'averages'};

      fits = load('data_fitting');
      params_labels = {'k_{AP}', '\alpha', 'k_{PA}', '\beta'};


      [junk, indxs] = ismember(targets, fits.all_data(:,1));
      all_data = fits.all_data(indxs,:);

      all_pts = cat(1, all_data{1:end-1,2});
      all_pts = [all_pts(:,2:5) 2*all_pts(:,end-3) all_pts(:,1)./all_pts(:,end)];
      %all_pts = [all_pts(:,1).*all_pts(:,4) all_pts(:,2) all_pts(:,3).*all_pts(:,2) all_pts(:,4:end)];

      outliers = pcout(all_pts(:, 1:4), [0.33 25], [0.25 1]);
      outliers_fraction = sum(outliers)/length(outliers);

      intercept = ones(sum(~outliers), 1);

      figure;
      subplot(2,2,1);
      myregress([intercept all_pts(~outliers,end)], all_pts(~outliers,1));
      subplot(2,2,2);
      myregress([intercept all_pts(~outliers,end)], all_pts(~outliers,2));
      subplot(2,2,3);
      myregress([intercept all_pts(~outliers,end)], all_pts(~outliers,3));
      subplot(2,2,4);
      myregress([intercept all_pts(~outliers,end)], all_pts(~outliers,4));
      mtit('Normalized score regress');

      figure;
      subplot(2,3,1);
      myregress([intercept all_pts(~outliers,5)], all_pts(~outliers,1));
      xlabel('Length of the embryo');
      ylabel(params_labels{1})
      subplot(2,3,2);
      myregress([intercept all_pts(~outliers,5)], all_pts(~outliers,2));
      xlabel('Length of the embryo');
      ylabel(params_labels{2})
      subplot(2,3,3);
      myregress([intercept all_pts(~outliers,5)], all_pts(~outliers,3));
      xlabel('Length of the embryo');
      ylabel(params_labels{3})
      subplot(2,3,4);
      myregress([intercept all_pts(~outliers,5)], all_pts(~outliers,4));
      xlabel('Length of the embryo');
      ylabel(params_labels{4})
      mtit('Embryo length regress');

      subplot(2,3,1)
      set(gca, 'YScale', 'log', 'YLim', [1e-5 100]);
      subplot(2,3,2)
      set(gca, 'YScale', 'log', 'YLim', [1e-1 10]);
      subplot(2,3,3)
      set(gca, 'YScale', 'log', 'YLim', [1e-5 100]);
      subplot(2,3,4)
      set(gca, 'YScale', 'log', 'YLim', [1e-1 10]);

      [m_all, s_all] = mymean(all_pts(:, 1:4));
      [m, s] = mymean(all_pts(~outliers, 1:4));

      disp(['All: ' num2str(m_all) ' +- ' num2str(s_all)]);
      disp(['Bests: ' num2str(m) ' +- ' num2str(s) ' (' num2str(1-outliers_fraction) ')']);

      all_pts = NaN(0, 6);

      for i=1:length(targets)-1
        pts = all_data{i,2};
        pts = [pts(:,2:5) 2*pts(:,end-3) pts(:,1)./pts(:,end)];

        %pts = [pts(:,1).*pts(:,4) pts(:,2) pts(:,3).*pts(:,2) pts(:,4:end)];
        goods = ~outliers(1:size(pts,1), 1);
        outliers = outliers(size(pts,1)+1:end, 1);

        [means, stds] = mymean(pts(goods,:));

        if (i==1)
          hfig(1) = figure('units','normalized','outerposition',[0 0 1 1])
        else
          figure(hfig(1));
        end
        subplot(2,3,1);hold on
        xlabel('Score');
        ylabel(params_labels{1})
        scatter(pts(goods,end), pts(goods,1), 'o', 'MarkerEdgeColor', colors(i,:));
        scatter(pts(~goods,end), pts(~goods,1), '+', 'MarkerEdgeColor', colors(i,:));
        errorbarxy(means(end), means(1), stds(end), stds(1), {colors(i+length(targets),:), colors(i+length(targets),:), colors(i+length(targets),:)});
        subplot(2,3,2);hold on
        xlabel('Score');
        ylabel(params_labels{2})
        scatter(pts(goods,end), pts(goods,2), 'o', 'MarkerEdgeColor', colors(i,:));
        scatter(pts(~goods,end), pts(~goods,2), '+', 'MarkerEdgeColor', colors(i,:));
        errorbarxy(means(end), means(2), stds(end), stds(2), {colors(i+length(targets),:), colors(i+length(targets),:), colors(i+length(targets),:)});
        subplot(2,3,3);hold on
        xlabel('Score');
        ylabel(params_labels{3})
        scatter(pts(goods,end), pts(goods,3), 'o', 'MarkerEdgeColor', colors(i,:));
        scatter(pts(~goods,end), pts(~goods,3), '+', 'MarkerEdgeColor', colors(i,:));
        errorbarxy(means(end), means(3), stds(end), stds(3), {colors(i+length(targets),:), colors(i+length(targets),:), colors(i+length(targets),:)});
        subplot(2,3,4);hold on
        xlabel('Score');
        ylabel(params_labels{4})
        scatter(pts(goods,end), pts(goods,4), 'o', 'MarkerEdgeColor', colors(i,:));
        scatter(pts(~goods,end), pts(~goods,4), '+', 'MarkerEdgeColor', colors(i,:));
        errorbarxy(means(end), means(4), stds(end), stds(4), {colors(i+length(targets),:), colors(i+length(targets),:), colors(i+length(targets),:)});
        if (i==length(targets)-1)
          subplot(2,3,1)
          set(gca, 'YScale', 'log', 'YLim', [1e-5 100]);
          drawnow
          set(gca, 'YTickLabel', roundn(get(gca,'ytick'), -5));
          subplot(2,3,2)
          set(gca, 'YScale', 'log', 'YLim', [1e-1 10]);
          drawnow
          set(gca, 'YTickLabel', roundn(get(gca,'ytick'), -4));
          subplot(2,3,3)
          set(gca, 'YScale', 'log', 'YLim', [1e-5 100]);
          drawnow
          set(gca, 'YTickLabel', roundn(get(gca,'ytick'), -5));
          subplot(2,3,4)
          set(gca, 'YScale', 'log', 'YLim', [1e-1 10]);
          drawnow
          set(gca, 'YTickLabel', roundn(get(gca,'ytick'), -4));
        end

        if (i==1)
          mtit(['Normalized score : ' num2str(outliers_fraction)])
          hfig(2) = figure('units','normalized','outerposition',[0 0 1 1])
        else
          figure(hfig(2));
        end
        subplot(2,3,1);hold on
        xlabel('Length of the embryo');
        ylabel(params_labels{1})
        scatter(pts(goods,5), pts(goods,1), 'o', 'MarkerEdgeColor', colors(i,:));
        scatter(pts(~goods,5), pts(~goods,1), '+', 'MarkerEdgeColor', colors(i,:));
        errorbarxy(means(5), means(1), stds(5), stds(1), {colors(i+length(targets),:), colors(i+length(targets),:), colors(i+length(targets),:)});
        subplot(2,3,2);hold on
        xlabel('Length of the embryo');
        ylabel(params_labels{2})
        scatter(pts(goods,5), pts(goods,2), 'o', 'MarkerEdgeColor', colors(i,:));
        scatter(pts(~goods,5), pts(~goods,2), '+', 'MarkerEdgeColor', colors(i,:));
        errorbarxy(means(5), means(2), stds(5), stds(2), {colors(i+length(targets),:), colors(i+length(targets),:), colors(i+length(targets),:)});
        subplot(2,3,3);hold on
        xlabel('Length of the embryo');
        ylabel(params_labels{3})
        scatter(pts(goods,5), pts(goods,3), 'o', 'MarkerEdgeColor', colors(i,:));
        scatter(pts(~goods,5), pts(~goods,3), '+', 'MarkerEdgeColor', colors(i,:));
        errorbarxy(means(5), means(3), stds(5), stds(3), {colors(i+length(targets),:), colors(i+length(targets),:), colors(i+length(targets),:)});
        subplot(2,3,4);hold on
        xlabel('Length of the embryo');
        ylabel(params_labels{4})
        scatter(pts(goods,5), pts(goods,4), 'o', 'MarkerEdgeColor', colors(i,:));
        scatter(pts(~goods,5), pts(~goods,4), '+', 'MarkerEdgeColor', colors(i,:));
        errorbarxy(means(5), means(4), stds(5), stds(4), {colors(i+length(targets),:), colors(i+length(targets),:), colors(i+length(targets),:)});
        if (i==length(targets)-1)
          subplot(2,3,1)
          set(gca, 'YScale', 'log', 'YLim', [1e-5 100]);
          drawnow
          set(gca, 'YTickLabel', roundn(get(gca,'ytick'), -5));
          subplot(2,3,2)
          set(gca, 'YScale', 'log', 'YLim', [1e-1 10]);
          drawnow
          set(gca, 'YTickLabel', roundn(get(gca,'ytick'), -4));
          subplot(2,3,3)
          set(gca, 'YScale', 'log', 'YLim', [1e-5 100]);
          drawnow
          set(gca, 'YTickLabel', roundn(get(gca,'ytick'), -5));
          subplot(2,3,4)
          set(gca, 'YScale', 'log', 'YLim', [1e-1 10]);
          drawnow
          set(gca, 'YTickLabel', roundn(get(gca,'ytick'), -4));
        end

        if (i==1)
          mtit('Egg length')
          hfig(3) = figure('units','normalized','outerposition',[0 0 1 1]);
        else
          figure(hfig(3));
        end

        count = 1;
        for k=1:3
          for l=k+1:4
            subplot(2,3,count);hold on
            xlabel(params_labels{k});
            ylabel(params_labels{l})
            scatter(pts(goods,k), pts(goods,l), 'o', 'MarkerEdgeColor', colors(i,:));
            %errorbarxy(means(k), means(l), stds(k), stds(l), {colors(i+length(targets),:), colors(i+length(targets),:), colors(i+length(targets),:)});
            count = count + 1;

            if (i==length(targets)-1)
              %yl = ylim;
              %xl = xlim;
              %ylim([0 yl(2)]);
              %xlim([0 xl(2)]);
              set(gca, 'XScale', 'log', 'YScale', 'log');

              if (k==1)
                set(gca, 'XLim', [1e-5 1], 'XTick', 10.^[-5:0])
              elseif (k==2)
                xl = xlim;
                set(gca, 'XLim', [0.5 10]);
              end
              if (l==2 || l==4)
                yl = ylim;
                set(gca, 'YLim', [0.5 10]);
              end

              drawnow
              set(gca, 'YTickLabel', roundn(get(gca,'ytick'), -6));
              set(gca, 'XTickLabel', roundn(get(gca,'xtick'), -6));
            end
          end
        end

        all_pts = [all_pts; [pts(goods,1:end-1) ones(sum(goods),1)*i]];

        if (i==1)
          mtit('Parameter correlation')
        end
      end

      %{
      intercept = ones(size(all_pts, 1), 1);
      figure;
      subplot(2,3,1);hold on;
      myregress([intercept all_pts(:,end-1)], all_pts(:,1));
      xlabel('Length of the embryo');
      ylabel(params_labels{1})
      subplot(2,3,2);hold on;
      myregress([intercept all_pts(:,end-1)], all_pts(:,2));
      xlabel('Length of the embryo');
      ylabel(params_labels{2})
      subplot(2,3,3);hold on;
      myregress([intercept all_pts(:,end-1)], all_pts(:,3));
      xlabel('Length of the embryo');
      ylabel(params_labels{3})
      subplot(2,3,4);hold on;
      myregress([intercept all_pts(:,end-1)], all_pts(:,4));
      xlabel('Length of the embryo');
      ylabel(params_labels{4})
      mtit('Egg length regress')
      %}

      figure;
      subplot(2,3,1);hold on;
      boxplot(all_pts(:,1), all_pts(:,end), 'label', all_data(1:end-1,1));
      ylabel(params_labels{1})
      set(gca, 'YScale', 'log');
      set(gca, 'YLim', [1e-5 1], 'YTick', 10.^[-5:0])
      set(gca, 'YTickLabel', roundn(get(gca,'ytick'), -6));
      subplot(2,3,2);hold on;
      boxplot(all_pts(:,2), all_pts(:,end), 'label', all_data(1:end-1,1));
      ylabel(params_labels{2})
      set(gca, 'YScale', 'log');
      set(gca, 'YLim', [0.5 10]);
      set(gca, 'YTickLabel', roundn(get(gca,'ytick'), -6));
      subplot(2,3,3);hold on;
      boxplot(all_pts(:,3), all_pts(:,end), 'label', all_data(1:end-1,1));
      ylabel(params_labels{3})
      set(gca, 'YScale', 'log');
      set(gca, 'YLim', [1e-3 1], 'YTick', 10.^[-3:0])
      set(gca, 'YTickLabel', roundn(get(gca,'ytick'), -6));
      subplot(2,3,4);hold on;
      boxplot(all_pts(:,4), all_pts(:,end), 'label', all_data(1:end-1,1));
      ylabel(params_labels{4})
      set(gca, 'YScale', 'log');
      set(gca, 'YLim', [0.5 10]);
      set(gca, 'YTickLabel', roundn(get(gca,'ytick'), -6));

      [H,P] = myttest(all_pts(:,1), all_pts(:,end))
      [H,P] = myttest(all_pts(:,2), all_pts(:,end))
      [H,P] = myttest(all_pts(:,3), all_pts(:,end))
      [H,P] = myttest(all_pts(:,4), all_pts(:,end))

      alpha = 0.05;
      bounds = prctile(all_pts, [alpha 0.5 1-alpha]*100);
      is_log = true(1, 4);

      figure;
      for i=1:4

        if (is_log(i))
          bound = max(max(abs(log10(bsxfun(@rdivide, bounds(:,logical(is_log)), bounds(2,logical(is_log)))))*1.75));
          pos = logspace(-bound, bound, 16)*bounds(2,i);
        else
          med = bounds(2,i);
          bound = max(abs(bounds([1 3],i) - med)*1.75);
          pos = linspace(med - bound, med + bound, 16);
        end

        %log_pos = logspace(-0.3, 0.3, 16);
        %tmp_pos = [-Inf log_pos(2:end-1) Inf];
        tmp_pos = [-Inf pos(2:end-1) Inf];

        nhist = histc(all_pts(:,i), tmp_pos);

        subplot(2,3,i);hold on
        %bar(log_pos, nhist(:,i),'histc');
        %bar(log_pos, nhist(:,i));
        bar(pos, nhist, 'histc');
        %if (i>1&&i<=5)
        if (is_log(i))
          set(gca,'XScale', 'log', 'XLim', pos([1 end]), 'XTick', diff(pos)/2+pos(1:end-1));
          set(gca, 'XTickLabel',roundn(get(gca,'XTick'), -5));
        else
          set(gca, 'XLim', pos([1 end]), 'XTick', diff(pos)/2+pos(1:end-1));
        end

        %hist(best_vals(~outliers,i), 16);
      end

      keyboard

    case 7.2
      vals = group_ml_results('LatestFits/adr-kymo-*_evol.dat', {'combine_data'}, {'type', '1056-temps-all'; 'fitting_type', 'sample'; 'fit_relative', true});

      for i = 1:size(vals, 1)
        best_indx = -1;
        date = 0;
        for j = 1:size(vals{i,2}, 1)
          curr_date = datenum(datevec(vals{i,2}{j,2}(end).time));
          if (date - curr_date < 0)
            date = curr_date;
            best_indx = j;
          end
        end

        if (best_indx > 0)
          vals{i,2} = vals{i,2}(best_indx,:);
        end
      end

      vals2 = extract_model_parameters(vals);
      sensitivity_analysis(vals2{1,2}{1,2}.evolution, 'best_model');

      [rC, C, rH, H] = correlation_matrix(vals2{2,2}{1,2}.evolution);

      figure;imagesc(rC, [-1 1]);
      colormap(flipud(redgreencmap));
      colorbar

      keyboard

    case 8.1
      colors = [254 217 166; ...
                251 180 174; ...
                222 203 228; ...
                255 127 0; ...
                228 26 28; ...
                152 78 163]/255;

      data = load('simul_optimized.mat');
      targets = {'good_ani2.txt';'good_24.txt';'good_c27d91.txt'};

      for p=1:length(data.all_simul)
        all_data = cell(3, 3);
        for f=1:size(data.all_simul{p}, 1)
          files = data.all_simul{p}{f,3};
          nfiles = size(files,1);

          all_maint = cell(nfiles, 1);
          centers = NaN(1, nfiles);
          all_params = NaN(nfiles, 6);

          for i=1:nfiles
            all_maint{i} = files{i,1}.';
            centers(i) = round((data.all_simul{p}{f,2}(i,1)/2)/opts_expansion.quantification.resolution);
            tmp_params = data.all_simul{p}{f,2}(i,1:5);
            tmp_params(1) = min(tmp_params(1), tmp_params(2)-1e-2);

            tmp_params = [tmp_params 2*tmp_params(3)];
            tmp_params(3) = ellipse_circum(tmp_params(3:5).', tmp_params(2), true);

            all_params(i,:) = tmp_params;
            disp([num2str(i) '/' num2str(nfiles)]);
          end
          [stack, shift] = stack_images(all_maint, centers, 0.5);
          center = centers(1) + shift(1);

          all_data{f,1} = stack;
          all_data{f,2} = all_params;
          all_data{f,3} = center;
        end

        [junk, indxs] = ismember(targets, data.all_simul{p}(:,1));
        all_data = all_data(indxs,:);
        pos = cell(1,3);

        figure('Name', num2str(p));
        for i=1:3
          stack = all_data{i,1};
          center = all_data{i,3};

          [m24, s24] = mymean(stack, 3);

          max_val = mean(m24(1:ntails));
          min_val = mean(m24(end-ntails+1:end));
          norm_factor = 1 / (max_val - min_val);

          m24 = (m24 - min_val) *norm_factor;
          s24 = s24 * norm_factor;
          stack = (stack - min_val) * norm_factor;

          all_data{i,1} = stack;

          pos{i} = ([1:size(stack,1)]-center)*opts_expansion.quantification.resolution;;

          subplot(2,3,i)
          hold on;
          plot(pos{i}, squeeze(stack),  'Color', colors(i, :));
          plot(pos{i}, m24+s24, 'Color', colors(i+3,:));
          plot(pos{i}, m24-s24, 'Color', colors(i+3,:));
          plot(pos{i}, m24, 'Color', colors(i+3,:), 'LineWidth', 2);
          ylim([0 1.2]);
          xlim([-20 20]);
          xlabel('Distance to the boundary (µm)');
          ylabel('GFP intensity (a.u.)')
          switch i
            case 1
              title('ani-2(RNAi)')
            case 2
              title('WT')
            case 3
              title('C27D9.1(RNAi)')
          end

          subplot(2,3,4);hold on
          plot(pos{i}, squeeze(stack),  'Color', colors(i, :));
          plot(pos{i}, m24, 'Color', colors(i+3, :), 'LineWidth', 2);
        end
        ylim([0 1.2]);
        xlim([-20 20]);

        [all_stacks, all_shift] = stack_images(all_data(:,1), cat(1,all_data{:,3}), 0.5);
        [m_all, s_all] = mymean(all_stacks, 3);
        pos_all = ([1:size(all_stacks,1)]-(all_data{1,3}+all_shift(1)))*opts_expansion.quantification.resolution;;

        subplot(2,3,4);hold on;
        %plot(pos_all, m_all, 'Color', [83 83 83]/255, 'LineWidth', 2);
        xlabel('Distance to the boundary (µm)');
        ylabel('GFP intensity (a.u.)')
        title('Average profiles')
        ylim([0 1.2]);
        xlim([-20 20]);

        subplot(2,3,5)
        hold on;
        plot(pos_all, squeeze(all_stacks), 'Color', [189 189 189]/255);
        plot(pos_all, m_all+s_all, 'Color', [83 83 83]/255);
        plot(pos_all, m_all-s_all, 'Color', [83 83 83]/255);
        plot(pos_all, m_all, 'Color', [83 83 83]/255, 'LineWidth', 2);
        xlabel('Distance to the boundary (µm)');
        ylabel('GFP intensity (a.u.)')
        title('Overall average profile')
        ylim([0 1.2]);
        xlim([-20 20]);

        disp([adaptive_neyman(squeeze(all_data{1,1}), squeeze(all_data{2,1})), ...
              adaptive_neyman(squeeze(all_data{1,1}), squeeze(all_data{3,1})), ...
              adaptive_neyman(squeeze(all_data{2,1}), squeeze(all_data{3,1}))]);

        avgs = NaN(3, 4);
        figure('Name', num2str(p));
        for i=1:3
          all_params = all_data{i,2};
          intercept = ones(size(all_params,1), 1);

          [means, stds] = mymean([all_params(:,6) all_params(:,1)./all_params(:,2)]);

          subplot(2,3,i);hold on
          title(['N=' num2str(length(intercept))])
          xlabel('Egg length (µm)');
          ylabel('Fraction (a.u.)')
          ylim([0 1]);
          xlim([30 80]);
          errorbarxy(means(1), means(2), stds(1), stds(2)); hold on
          myregress([intercept all_params(:,6)], all_params(:,1)./all_params(:,2), colors(i, :), 'dx');

          subplot(2,3,4);hold on
          ylim([0 1]);
          xlim([30 80]);
          errorbarxy(means(1), means(2), stds(1), stds(2)); hold on
          myregress([intercept all_params(:,6)], all_params(:,1)./all_params(:,2), colors(i, :), 'dx');
        end
        title('ani-2, wt, C27D9.1')
        xlabel('Egg length (µm)');
        ylabel('Fraction (a.u.)')
        ylim([0 1]);
        xlim([30 80]);

        all_params = cat(1, all_data{:,2});
        intercept = ones(size(all_params,1), 1);
        [means, stds] = mymean([all_params(:,6) all_params(:,1)./all_params(:,2)]);

        subplot(2,3,5);hold on
        ylim([0 1]);
        xlim([30 80]);
        errorbarxy(means(1), means(2), stds(1), stds(2));hold on
        myregress([intercept all_params(:,6)], all_params(:,1)./all_params(:,2), black, 'dx');
        title('All together')
        xlabel('Egg length (µm)');
        ylabel('Fraction (a.u.)')

        mtit('Simulation of membrane length fraction')

        [ratios, surfaces, volumes] = surface2volume(all_params(:,3:5).');
        heights = find_cap_size(all_params(:,3:5).', all_params(:,1).');
        [caps_surface, caps_volume] = spheroidal_cap(all_params(:,3:5).', heights);

        figure('Name', num2str(p));

        subplot(2,3,1);
        ylim([0 1]);
        xlim([30 80]);
        myregress([intercept all_params(:,6)], (all_params(:,3)-heights.')./(2*all_params(:,3)), black, 'dx');
        title('Major axis fraction')
        xlabel('Egg length (µm)');
        ylabel('Fraction (a.u.)')

        subplot(2,3,2);
        ylim([0 1]);
        xlim([30 80]);
        myregress([intercept all_params(:,6)], caps_surface./surfaces, black, 'dx');
        title('Membrane surface fraction')
        xlabel('Egg length (µm)');
        ylabel('Fraction (a.u.)')

        subplot(2,3,3);
        ylim([0 1]);
        xlim([30 80]);
        myregress([intercept all_params(:,6)], caps_volume./volumes, black, 'dx');
        title('Cell volume fraction')
        xlabel('Egg length (µm)');
        ylabel('Fraction (a.u.)')

        mtit('Alternative fractions')

        subplot(2,3,4);
        myregress([intercept all_params(:,6)], all_params(:,1), black, 'dx');
        title('Domain size vs egg length')
        xlabel('Egg length (\mu m)');
        ylabel('Domain length (\mu m)')

        subplot(2,3,5);
        myregress([intercept all_params(:,2)], all_params(:,1), black, 'dx');
        title('Domain size vs membrane length')
        xlabel('Membrane length (µm)');
        ylabel('Domain length (\mu m)')

        subplot(2,3,6);
        myregress(all_params(:,2), all_params(:,1), black, 'dx');
        title('Domain size vs membrane length, no intercept')
        xlabel('Membrane length (µm)');
        ylabel('Domain length (\mu m)')
        ylim([30 90])
        xlim([80 180])


        mtit('Simulations')
      end

      %keyboard

    case 8.2

      colors = [179 205 227; ...
                204 235 197; ...
                251 180 174; ...
                55 126 184; ...
                77 175 74; ...
                228 26 28]/255;

      data = load('simul_optimized.mat');
      targets = {'good_13.txt';'good_20.txt';'good_24.txt'};

      for p=1:length(data.all_simul)
        all_data = cell(3, 3);
        for f=1:size(data.all_simul{p}, 1)
          files = data.all_simul{p}{f,3};
          nfiles = size(files,1);

          all_maint = cell(nfiles, 1);
          centers = NaN(1, nfiles);
          all_params = NaN(nfiles, 6);

          for i=1:nfiles
            all_maint{i} = files{i,1}.';
            centers(i) = round((data.all_simul{p}{f,2}(i,1)/2)/opts_expansion.quantification.resolution);

            tmp_indx = NaN(size(thresh));
            fraction = files{i,2};
            for t=1:length(thresh)
              tmp_indx(t) = find(fraction>=thresh(t), 1, 'first');
            end

            all_params(i,:) = [data.all_simul{p}{f,2}(i,1:5) (tmp_indx(2)-tmp_indx(1))*10];
            disp([num2str(i) '/' num2str(nfiles)]);
          end
          [stack, shift] = stack_images(all_maint, centers, 0.5);
          center = centers(1) + shift(1);

          all_data{f,1} = stack;
          all_data{f,2} = all_params;
          all_data{f,3} = center;
        end

        [junk, indxs] = ismember(targets, data.all_simul{p}(:,1));
        temperatures = data.temperatures(indxs)
        all_data = all_data(indxs,:);
        pos = cell(1,3);

        all_times = NaN(0,2);
        avgs = NaN(1,3);

        figure('Name', num2str(p));hold on;
        for i=1:3
          boxplot(all_data{i,2}(:,6), 'colors',colors(i, :), 'position', temperatures(i), 'width', 2, 'jitter', 1);
          all_times = [all_times; [all_data{i,2}(:,6) ones(size(all_data{i,2},1),1)*i]];

          hs = findobj(gcf,'tag','Outliers');
          xc = get(hs,'XData');
          if (i>1)
            xc = xc{1};
          end

          if (isnan(xc))
            goods = true(size(all_data{i,2}, 1), 1);
          else
            yc = get(hs,'YData')

            if (i>1)
              yc = yc{1};
            end

            goods = ~(ismember(all_data{i,2}(:,6), yc));
          end

          avgs(i) = mymean(all_data{i,2}(goods,6));
        end
        %avgs = mymean(all_times(:,1), 1, all_times(:,2));
        plot(temperatures, avgs, '-ok');

        title('Polarity establishment')
        ylim([0 max(all_times(:,1)+10)]);
        xlim([10 26]);
        set(gca, 'XTick', [10:26], 'XTickLabel', [10:26]);
        [H,P] = myttest(all_times(:,1), all_times(:,end))

        figure('Name', num2str(p));
        for i=1:3
          stack = all_data{i,1};
          center = all_data{i,3};

          [m24, s24] = mymean(stack, 3);
          max_val = mean(m24(1:ntails));
          min_val = mean(m24(end-ntails+1:end));
          norm_factor = 1 / (max_val - min_val);

          m24 = (m24 - min_val) *norm_factor;
          s24 = s24 * norm_factor;
          stack = (stack - min_val) * norm_factor;

          all_data{i,1} = stack;

          pos{i} = ([1:size(stack,1)]-center)*opts_expansion.quantification.resolution;;

          subplot(2,3,i)
          hold on;
          plot(pos{i}, squeeze(stack),  'Color', colors(i, :));
          plot(pos{i}, m24+s24, 'Color', colors(i+3,:));
          plot(pos{i}, m24-s24, 'Color', colors(i+3,:));
          plot(pos{i}, m24, 'Color', colors(i+3,:), 'LineWidth', 2);
          ylim([0 1.2]);
          xlim([-20 20]);
          xlabel('Distance to the boundary (µm)');
          ylabel('GFP intensity (a.u.)')
          switch i
            case 1
              title('13°C')
            case 2
              title('20°C')
            case 3
              title('24°C')
          end

          subplot(2,3,4);hold on
          plot(pos{i}, squeeze(stack),  'Color', colors(i, :));
          %plot(pos{i}, m24+s24,  'Color', colors(i, :));
          %plot(pos{i}, m24-s24,  'Color', colors(i, :));
          plot(pos{i}, m24, 'Color', colors(i+3, :), 'LineWidth', 2);
        end
        ylim([0 1.2]);
        xlim([-20 20]);

        [all_stacks, all_shift] = stack_images(all_data(:,1), cat(1,all_data{:,3}), 0.5);
        [m_all, s_all] = mymean(all_stacks, 3);
        pos_all = ([1:size(all_stacks,1)]-(all_data{1,3}+all_shift(1)))*opts_expansion.quantification.resolution;;

        subplot(2,3,4);hold on;
  %      plot(pos_all, m_all+s_all, 'Color', [83 83 83]/255);
  %      plot(pos_all, m_all-s_all, 'Color', [83 83 83]/255);
        plot(pos_all, m_all, 'Color', [83 83 83]/255, 'LineWidth', 2);
        xlabel('Distance to the boundary (µm)');
        ylabel('GFP intensity (a.u.)')
        title('Average profiles')
        ylim([0 1.2]);
        xlim([-20 20]);

        subplot(2,3,5)
        hold on;
        plot(pos_all, squeeze(all_stacks), 'Color', [189 189 189]/255);
        plot(pos_all, m_all+s_all, 'Color', [83 83 83]/255);
        plot(pos_all, m_all-s_all, 'Color', [83 83 83]/255);
        plot(pos_all, m_all, 'Color', [83 83 83]/255, 'LineWidth', 2);
        xlabel('Distance to the boundary (µm)');
        ylabel('GFP intensity (a.u.)')
        title('Overall average profile')
        ylim([0 1.2]);
        xlim([-20 20]);

        disp([adaptive_neyman(squeeze(all_data{1,1}), squeeze(all_data{2,1})), ...
              adaptive_neyman(squeeze(all_data{1,1}), squeeze(all_data{3,1})), ...
              adaptive_neyman(squeeze(all_data{2,1}), squeeze(all_data{3,1}))]);
        
        figure('Name', num2str(p));
        for i=1:3
          all_params = all_data{i,2};
          intercept = ones(size(all_params,1), 1);
          [means, stds] = mymean([2*all_params(:,3) all_params(:,1)./all_params(:,2)]);

          subplot(2,3,i);hold on
          ylim([0 1]);
          xlim([30 80]);
          errorbarxy(means(1), means(2), stds(1), stds(2));hold on;
          myregress([intercept 2*all_params(:,3)], all_params(:,1)./all_params(:,2), colors(i, :));
          title(['N=' num2str(length(intercept))])
          xlabel('Egg length (µm)');
          ylabel('Fraction (a.u.)')

          subplot(2,3,4);hold on
          ylim([0 1]);
          xlim([30 80]);
          errorbarxy(means(1), means(2), stds(1), stds(2)); hold on
          myregress([intercept 2*all_params(:,3)], all_params(:,1)./all_params(:,2), colors(i, :));
        end
        title('13, 20, 24')
        xlabel('Egg length (µm)');
        ylabel('Fraction (a.u.)')
        ylim([0 1]);
        xlim([30 80]);

        all_params = cat(1, all_data{:,2});
        intercept = ones(size(all_params,1), 1);

        subplot(2,3,5);hold on
        ylim([0 1]);
        xlim([30 80]);
        myregress([intercept 2*all_params(:,3)], all_params(:,1)./all_params(:,2));
        title('All together')
        xlabel('Egg length (µm)');
        ylabel('Fraction (a.u.)')

        mtit('Membrane length fraction')
      end

%      keyboard

    case 9

      load('data_expansion');
      tmp_data = cat(1, all_data{:,3});
      all_names = tmp_data(:, 3);

      tmp_opts = get_struct('modeling');
      opts_goehring = load_parameters(tmp_opts, 'goehring.txt');
      opts_flow = load_parameters(opts_goehring, 'custom_flow.txt');

      orig_params = [opts_flow.diffusion_params; opts_flow.reaction_params];
      orig_params([4 12]) = orig_params([4 12]) ./ orig_params([13 5]);
      orig_params = orig_params(:).';

      all_scores = NaN(length(all_names), 5);
      all_offsets = all_scores;
      all_params = all_scores(:,1:3);

      vals = group_ml_results('LatestFits/adr-kymo-*_evol.dat', {'type';'simulation_parameters';'flow_size';'scale_flow';'fit_flow'}, {'parameter_set', 0});
      vals = extract_model_parameters(vals, true);

      nfits = size(vals,1);

      for i=1:nfits
        curr_val = vals{i,1}{1};

        indx = find(ismember(all_names, [curr_val.type '.mat']), 1);
        if (isempty(indx))
          continue;
        end

        if (~curr_val.scale_flow)
          if (curr_val.flow_size(2) == size(opts_goehring.advection_params, 2))
            sub_indx = 1;
          elseif (numel(curr_val.simulation_parameters)>4 && all(curr_val.simulation_parameters(4:5) == orig_params(4:5)))
            sub_indx = 2;
          else
            sub_indx = 3;
          end
        else
          if (~curr_val.fit_flow)
            sub_indx = 4;
            
            all_params(indx, 1) = vals{i,2}{1,2}.params.flow_scaling;
          else
            sub_indx = 5;

            all_params(indx, 2) = vals{i,2}{1,2}.params.flow_scaling;
            all_params(indx, 3) = vals{i,2}{1,2}.params.flow;
          end
        end

        all_scores(indx, sub_indx) = vals{i,2}{1,2}.score;
        all_offsets(indx, sub_indx) = vals{i,2}{1,2}.params.offset;
      end

      vals = group_ml_results('LatestFits/adr-kymo-*_evol.dat', {'parameter_set';'simulation_parameters';'flow_size';'scale_flow';'fit_flow'}, {'type', '1056-all-all'; 'scale_each_egg', true; 'extrapol_z', true});
      vals = extract_model_parameters(vals, true);

      nfits = size(vals, 1);
      all_scores2 = NaN(nfits, 5);
      all_params2 = cell(nfits, 2);

      for i=1:nfits
        curr_val = vals{i,1}{1};
        best = Inf;

        for j=1:size(vals{i,2}, 1)
          value = vals{i,2}{j,2};
          if (value.score < best)
            all_scores2(i, :) = [value.score curr_val.parameter_set curr_val.fit_flow curr_val.scale_flow curr_val.flow_size(2)];
            all_params2{i,1} = [value.params.rate value.params.flow value.params.flow_scaling];
            all_params2{i,2} = curr_val.simulation_parameters;

            best = value.score;
          end
        end
      end

      [junk, indx] = sort(all_scores2(:,1), 'descend');
      all_scores2 = all_scores2(indx, :);
      all_params2 = all_params2(indx, :);

      keyboard

  end

  return;
end

function [means, center, dwidth, full_width, f] = get_profile(domain, nframes, opts)

  opts_expansion = load_parameters(get_struct('ASSET'), 'domain_expansion.txt');
  prct_thresh = 5;
  has_noise = true;

  if (nargin == 3)
    mymovie = domain;
    [f, frac_width, full_width, domain] = domain_expansion(mymovie, opts);
    [h, w, nplanes] = size(domain);
    width = (w-1)/2;
  else
    has_noise = false;

    [h, w, nplanes] = size(domain);
    width = (w-1)/2;

    avg_domain = mymean(domain,3);

    bads = all(isnan(avg_domain),2);
    first = find(~bads, 1, 'first');
    last = find(~bads, 1, 'last');

    domain = domain(first:last,:,:);
    avg_domain = avg_domain(first:last,:);

    [f, frac_width, full_width] = domain_expansion(avg_domain, width + 1, last-first+1, opts_expansion);
  end

  domain = permute([domain(:, 1:width+1, :), domain(:,end:-1:width+1, :)], [2 1 3]);

  path = f*frac_width/opts_expansion.quantification.resolution;

  pos_mat = repmat([1:width].', 1, length(path));
  mask = bsxfun(@le, pos_mat, path.');
  mask = repmat(flipud(mask), [2, 1]);

  if (has_noise)
    noise = estimate_noise(domain);
    domain = min_max_domain(domain, path, 3*noise(:,2));
  end

  for i=1:nplanes
    img = domain(:,:,i);
    min_val = prctile(img(~mask), prct_thresh);
    max_val = prctile(img(mask), 100-prct_thresh);
    if (min_val > max_val)
      domain(:,:,i) = (img - max_val) / (min_val - max_val);
    else
      domain(:,:,i) = (img - min_val) / (max_val - min_val);
    end
  end

  domain = permute([domain(1:width+1, :, :); domain(end-1:-1:width+2, :, :)], [2 1 3]); 

  data = domain(end-nframes+1:end, :, :);
  data = [data(:,(width + 1):end,:); data(:,(width + 1):-1:1,:)];
  data = reshape(permute(data, [1 3 2]), [], width+1);
  [means, stds] = mymean(data, 1);

  dwidth = mymean(path(end-nframes+1:end));
  center = round(dwidth);
  dwidth = dwidth*opts_expansion.quantification.resolution;

  %figure;plot(means);
  %hold on;
  %scatter(center, means(center), 20, 'r')
  %ylim([0 1])

  if (nargout == 2)
    stds = center;
  end

  return;
end
