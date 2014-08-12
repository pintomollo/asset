function figs_msb(num)

  if (nargin == 0)
    num = 1;
  end
  opts_expansion = load_parameters(get_struct('ASSET'), 'domain_expansion.txt');
  colors = 'gbrcmyk';
  black = ones(1,3)*37/255;

  rescale_size = [400 700];
  thresh = [0 0.775];
  nframes = 25;
  ntails = 20;

  switch num
    case -1

      all_files = {'good_24.txt'; 'good_20.txt'; 'good_13.txt'; 'good_c27d91.txt'; 'good_ani2.txt'};
      [time, names] = get_manual_timing();
      all_data = NaN(0,4);

      for f = 1:length(all_files)
        files = textread(all_files{f}, '%s');
        nfiles = length(files);

        all_sizes = NaN(nfiles, 4);
        for i=1:nfiles
          load(files{i});
          good_time = ismember(names, files{i});
          pnm = time(good_time, 2);

          if (~isnan(pnm))
            cortex = mymovie.dic.cortex(pnm).carth;
            cortexr = realign(cortex, [1 1], mymovie.dic.centers(:,pnm), mymovie.dic.orientations(1,pnm));
            conv = convhull(cortexr(:,1),cortexr(:,2));
            [junk,axes_length,junk]=fit_ellipse(cortexr(conv,:));

            all_sizes(i,:) = [mymovie.dic.axes_length(:,pnm).' axes_length.'];
          end
          disp([num2str(i) '/' num2str(nfiles)]);
        end
        all_data = [all_data; all_sizes];
      end

      % Avg difference with egg size:
      % 1 + [3.4 1.6]/100

      keyboard

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

      load('data_expansion');

      convert = get_struct('z-correlation');

      all_params = cat(1, all_data{:,2});
      all_params(:, 5) = convert.bkg + convert.long_axis*all_params(:, 3) + convert.short_axis*all_params(:, 4);
      all_params(:,3) = ellipse_circum(all_params(:,3:5), all_params(:,2), true);

      all_sizes = all_params(:,3:5);
      all_sizes = [all_sizes all_params(:,2)*0.5 surface2volume(all_sizes.').'];

      data = load('data_fitting');
      all_params = cat(1, data.all_data{1:5,2});
      all_params = all_params(:, 2:5);

      opts = get_struct('modeling');
      opts = load_parameters(opts, 'goehring.txt');
      opts = load_parameters(opts, 'custom_flow.txt');

      flow = opts.advection_params;
      if (size(flow, 1) ~= opts.nparticles)
        [X, Y] = meshgrid([1:size(flow, 2)], 1+([0:opts.nparticles-1]*(size(flow, 1)-1)/(opts.nparticles-1)).');
        flow = bilinear_mex(flow, X, Y, [2 2]);
      end

      all_simul = cell(size(all_params, 1), 1);

      for p = 1:size(all_params, 1)
        opts.reaction_params([3 4 10 11]) = all_params(p, :);
        counter = 1;
        for f = 1:size(all_data, 1)
          for i=1:size(all_data{f,2}, 1)
            opts.axes_length = all_sizes(counter, 1:3).';

            opts.reaction_params(end, :) = all_sizes(counter, 4);
            opts.reaction_params(end-1,:) = all_sizes(counter, 5);

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

            disp([num2str(counter) '/' num2str(size(all_sizes,1)) '/' num2str(p)]);
            counter = counter + 1;
          end
        end

        all_simul{p} = all_data;
      end

      save('simul_expansion', 'all_simul', 'temperatures');

      keyboard

    case 0.3
      load('data_expansion')
      data = load('1056-all-all.mat');

      models = {'goehring', 'custom_flow','average_model', 'extended_model', 'full_model', 'final_model'};
      %models = {'test_model1', 'test_model2', 'final_model'};
      %models = {'test_model4', 'final_model'};
      all_simul = cell(length(models), 1);
      scores = NaN(length(models), 1);

      %[scores, results] = check_parameters('1056-all-all.mat', 'config_fitting', 'fit_kymo', 'config_modeling', 'full_model', 'start_with_best', false, 'aligning_type', 'best', 'init_noise', false);

      for p=1:length(models)
        [scores(p), results] = check_parameters('1056-all-all.mat', 'config_fitting', 'fit_kymo', 'config_modeling', models{p}, 'start_with_best', false, 'aligning_type', 'best', 'init_noise', false);

        count = 1;
        for f = 1:size(all_data, 1)
          for i = 1:size(all_data{f,2}, 1)

            domain = results{count, 1};

            domain = domain((end/2)+1:end, :).';
            domain = [domain domain(:,end-1:-1:1)];

            [profile, center, max_width, cell_width, path] = get_profile(domain, nframes);

            norig = length(all_data{f,3}{i,1})-1;
            nprofile = length(profile)-1;
            profile = interp1([0:nprofile], profile, [0:norig]*nprofile/norig);
            dwidth = mymean(results{count,2}(end-nframes+1:end));

            all_data{f,2}(i,1) = min(2*dwidth, all_data{f,2}(i,2));
            all_data{f,3}{i,1} = profile;
            all_data{f,3}{i,2} = [path(:) results{count, 2}(:)];

            disp([num2str(i) '/' num2str(size(all_data{f,2},1))]);

            count = count + 1;
          end
        end

        all_simul{p} = all_data;
      end
      temperatures = data.temperatures;

      save('simul_optimized', 'all_simul', 'temperatures', 'scores');

      keyboard

    case 0.4
      opts = get_struct('modeling');
      opts = load_parameters(opts, 'goehring.txt');

      %ratios = [1.56 [1.5:-0.1:1.3]];
      %ratios = [1.56:-0.1:1.2];
      ratios = [0.5:0.25:1.5];

      data = load('1056-24-all.mat');
      domains = cell(length(data), 1);

      flow = opts.advection_params;
      if (size(flow, 1) ~= opts.nparticles)
        [X, Y] = meshgrid([1:size(flow, 2)], 1+([0:opts.nparticles-1]*(size(flow, 1)-1)/(opts.nparticles-1)).');
        flow = bilinear_mex(flow, X, Y, [2 2]);
      end

      egg_size = opts.axes_length;
      cell_width = 2*max(data.pos);

      egg_size = sort(egg_size, 'descend');
      convert = get_struct('z-correlation');
      egg_size(3,:) = convert.bkg + convert.long_axis*egg_size(1,:) + convert.short_axis*egg_size(2,:);

      egg_size(1) = ellipse_circum(egg_size, cell_width, true);

      opts.axes_length = egg_size;

      opts.reaction_params(end-1,:) = surface2volume(opts.axes_length);
      opts.reaction_params(end, :) = 0.5*ellipse_circum(opts.axes_length);

      opts.boundaries = [0 opts.reaction_params(end,1)];
      opts.x_step = diff(opts.boundaries)/(opts.nparticles-1);

      corr_offset = 8;
      half = ((length(data.pos)-1)/2)+1;

      ml_params = [opts.diffusion_params; ...
                    opts.reaction_params];

      orig_params = ml_params;

      hfig = figure;
      for j=1:numel(ml_params)-3
        for i=1:length(ratios)
          ml_params = orig_params;
          ml_params(j) = ml_params(j)*ratios(i);

          opts.diffusion_params = ml_params(1,:);
          opts.reaction_params = ml_params(2:end,:);
          %opts.reaction_params(5,1) = ratios(i);

          x0 = opts.init_func(opts, false);


          [res, t_pos] = simulate_model_mix(x0, ml_params, opts.x_step, opts.tmax*0.75, opts.time_step, opts.output_rate, flow, opts.user_data, opts.max_iter);

          %keyboard

          %res = res((end/2)+1:end, :);
          %res = flipud(interp1q([0:size(res, 1)].'*opts.x_step, flipud(res), data.pos(half:end).'));

          res = flipud(interp2(flipud(res), [1:size(data.ground_truth, 1)]+corr_offset, data.pos(half:end).'/opts.x_step + 1));

          %domain = domain((end/2)+1:end, :).';
          domain = res.';
          domain = imnorm([domain domain(:,end-1:-1:1)]);

          domains{i} = domain;

          %figure;imagesc(domain);title(num2str(ratios(i)));
          imagesc(domain);title(num2str([j ratios(i)]));
          colormap(blueredmap)
          plot2svg(['PNG/enrichment_sensitivity_' num2str(j) '_' num2str(i) '.svg'], hfig);

          %figure;
          %find_kymograph('1056-24-all.mat', opts, 'config_fitting', 'fit_kymo', 'display', true, 'aligning_type', 'best', 'init_noise', 0)

          disp([num2str(i) '/' num2str(length(ratios))]);
        end
        disp(['---' num2str(j) '/' num2str(numel(ml_params)+1)]);
      end

      j=j+1;
      for i=1:length(ratios)
        ml_params = orig_params;

        opts.diffusion_params = ml_params(1,:);
        opts.reaction_params = ml_params(2:end,:);
        %opts.reaction_params(5,1) = ratios(i);

        x0 = opts.init_func(opts, false);


        [res, t_pos] = simulate_model_mix(x0, ml_params, opts.x_step, opts.tmax*0.75, opts.time_step, opts.output_rate, flow*ratios(i), opts.user_data, opts.max_iter);

        %keyboard

        %res = res((end/2)+1:end, :);
        %res = flipud(interp1q([0:size(res, 1)].'*opts.x_step, flipud(res), data.pos(half:end).'));

        res = flipud(interp2(flipud(res), [1:size(data.ground_truth, 1)]+corr_offset, data.pos(half:end).'/opts.x_step + 1));

        %domain = domain((end/2)+1:end, :).';
        domain = res.';
        domain = imnorm([domain domain(:,end-1:-1:1)]);

        domains{i} = domain;

        %figure;imagesc(domain);title(num2str(ratios(i)));
        imagesc(domain);title(num2str([j ratios(i)]));
        colormap(blueredmap)
        plot2svg(['PNG/enrichment_sensitivity_' num2str(j) '_' num2str(i) '.svg'], hfig);

        %figure;
        %find_kymograph('1056-24-all.mat', opts, 'config_fitting', 'fit_kymo', 'display', true, 'aligning_type', 'best', 'init_noise', 0)

        disp([num2str(i) '/' num2str(length(ratios))]);
      end
      disp(['---' num2str(j) '/' num2str(numel(ml_params)+1)]);

      %[profile, center, max_width, cell_width, path] = get_profile(domain, nframes);
      %norig = length(all_data{f,3}{i,1})-1;
      %nprofile = length(profile)-1;
      %profile = interp1([0:nprofile], profile, [0:norig]*nprofile/norig);

      %all_data{f,2}(i,1) = 2*max_width * opts.x_step/opts_expansion.quantification.resolution;
      %all_data{f,3}{i,1} = profile;
      %all_data{f,3}{i,2} = path;

      keyboard

    case 0.5
      load('1056-24-040311_1_');
      nframes = size_data(mymovie.dic);

      for nimg=1:nframes
        img = imnorm(double(load_data(mymovie.dic,nimg)));
        img = mask_neighbors(img, mymovie.dic.centers(:,nimg), mymovie.dic.axes_length(:,nimg), mymovie.dic.orientations(1,nimg), mymovie.dic.neighbors(nimg), opts);
        img = realign(img,rescale_size,mymovie.dic.centers(:,nimg),mymovie.dic.orientations(1,nimg));
        imwrite(img, ['PNG/Movie/DIC-' num2str(nimg) '.tif'], 'tif');

        img = imnorm(double(load_data(mymovie.cortex,nimg)));
        img = mask_neighbors(img, mymovie.markers.centers(:,nimg), mymovie.markers.axes_length(:,nimg), mymovie.markers.orientations(1,nimg), mymovie.markers.neighbors(nimg), opts);
        img = realign(img,rescale_size,mymovie.markers.centers(:,nimg),mymovie.markers.orientations(1,nimg));
        imwrite(img, ['PNG/Movie/PH-' num2str(nimg) '.tif'], 'tif');

        img = imnorm(double(load_data(mymovie.data,nimg)));
        img = mask_neighbors(img, mymovie.data.centers(:,nimg), mymovie.data.axes_length(:,nimg), mymovie.data.orientations(1,nimg), mymovie.data.neighbors(nimg), opts);
        img = realign(img,rescale_size,mymovie.data.centers(:,nimg),mymovie.data.orientations(1,nimg));
        imwrite(img, ['PNG/Movie/PAR2-' num2str(nimg) '.tif'], 'tif');

        disp([num2str(nimg) '/' num2str(nframes)]);
      end


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

      %%%%%%%%%%%%%%%%%%%%%%%% Importing from ImageJ the other data

      files = dir('signal_quantif/*.tif');
      nfiles = length(files);
      names = cell(nfiles, 1);

      prct = 50;
      all_vals = NaN(0, 5);
      all_noises = NaN(0,5);

      pos = NaN(1,4);
      ngroups = NaN(1,4);
      for i=1:nfiles
        fname = ['signal_quantif/' files(i).name];
        names{i,1} = files(i).name(1:find(files(i).name=='_')-1);

        [nimg, size_img] = size_data(fname);
        roi_name = dir(['signal_quantif/' names{i,1} '*.zip']);
        rois = ReadImageJROI(fullfile(pwd, ['signal_quantif/' roi_name(1).name]));

        pixels = cell(nimg, 1);
        noises = NaN(nimg, 4);
        data = NaN(nimg,4);

        for n=1:nimg
          img = double(load_data(fname, n));
          img(img == 0) = NaN;
          noises(n,:) = estimate_noise(img);
          indx = (n-1)*4;

          for j = 1:4
            curr_roi = (rois{indx+j}.vfShapes);
            mask = false(size_img);

            for k = 1:length(curr_roi)
              mask = mask | roipoly(mask, curr_roi{k}(:,1), curr_roi{k}(:,2));
            end

            pts = img(mask);
            thresh = prctile(pts, prct);
            data(n, j) = mean(pts(pts > thresh));
            pos(j) = rois{indx+j}.vnRectBounds(2);
            ngroups(j) = k;
          end

          cyto = pos(ngroups >= 5);
          pts_indx = find(ngroups>=5);
          [junk, tmp_indx] = max(cyto);
          %if (cyto(1) < cyto(2))
          %  pts_indx = pts_indx([2 1]);
          %end
          pts_indx = [find(ngroups < 5) 1 pts_indx(tmp_indx) 1];

          %cortex = pos(ngroups < 5);
          %if (cortex(1) < cortex(2))
          %  pts_indx = [fliplr(find(ngroups < 5)) pts_indx];
          %else
          %  pts_indx = [find(ngroups < 5) pts_indx];
          %end

          data(n, :) = data(n, pts_indx);

          disp([num2str(n) '/' num2str(nimg)]);
        end

        all_vals = [all_vals; [data ones(nimg, 1)*i]];
        all_noises = [all_noises; [noises ones(nimg, 1)*i]];
        disp([num2str(i) '/' num2str(nfiles)]);
      end

      noises_grp = mymean(all_noises, 1, all_noises(:,end));
      noise = noises_grp(all_vals(:,end),:);

      val = (all_vals(:,1) - all_vals(:,3)) ./ noise(:,2);
      m = mymean(val, 1, all_vals(:,end));
      val = val / m(1);
      figure;boxplot(val, all_vals(:,end))
      title('SNR GFP')

      [H,p] = myttest(val, all_vals(:,end))

      %%%%%%%%%%%%%% Direct import from ImageJ to get percentiles

      files = dir('stainings/*.tif');
      nfiles = length(files);
      names = cell(nfiles, 1);

      prct = 50;
      all_vals = NaN(0, 5);
      all_noises = NaN(0,5);
      %all_vals2 = NaN(0, 5);
      %all_data = cell(nfiles, 2);

      file = dir('stainings/1056*.tif');

      fname = ['stainings/' file(1).name];
      roi_name = dir(['stainings/1056*.zip']);
      rois = ReadImageJROI(fullfile(pwd, ['stainings/' roi_name(1).name]));

      img = double(load_data(fname, 1));
      img(img == 0) = NaN;
      noise = estimate_noise(img);

      ranges = NaN(1, 4);
      for j=1:4
        ranges = [min(ranges(1:2), rois{j}.vnRectBounds(1:2)) max(ranges(3:4), rois{j}.vnRectBounds(3:4))];
      end
      cntr = round(mean([ranges(1:2);ranges(3:4)]));
      cntr = cntr([2 1]);

      figure;
      imshow(realign(imnorm(img-noise(1)),rescale_size,cntr,0));
      hold on

      for j = 1:4
        curr_roi = (rois{j}.vfShapes);

        for k = 1:length(curr_roi)
          plot(curr_roi{k}(:,1)-cntr(1)+rescale_size(2)/2, curr_roi{k}(:,2)-cntr(2)+rescale_size(1)/2);
        end
      end

      pos = NaN(1,4);
      ngroups = NaN(1,4);
      for i=1:nfiles
        fname = ['stainings/' files(i).name];
        names{i,1} = files(i).name(1:find(files(i).name=='_')-1);

        [nimg, size_img] = size_data(fname);
        roi_name = dir(['stainings/' names{i,1} '*.zip']);
        rois = ReadImageJROI(fullfile(pwd, ['stainings/' roi_name(1).name]));

        %all_pts = NaN(nimg, 3);

        pixels = cell(nimg, 1);
        noises = NaN(nimg, 4);
        data = NaN(nimg,4);

        for n=1:nimg
          img = double(load_data(fname, n));
          img(img == 0) = NaN;
          noises(n,:) = estimate_noise(img);
          indx = (n-1)*4;

          for j = 1:4
            curr_roi = (rois{indx+j}.vfShapes);
            mask = false(size_img);

            for k = 1:length(curr_roi)
              mask = mask | roipoly(mask, curr_roi{k}(:,1), curr_roi{k}(:,2));
            end

            pts = img(mask);
            thresh = prctile(pts, prct);
            data(n, j) = mean(pts(pts > thresh));
            pos(j) = rois{indx+j}.vnRectBounds(2);
            ngroups(j) = k;
          end

          cyto = pos(ngroups >= 5);
          pts_indx = find(ngroups>=5);
          if (cyto(1) < cyto(2))
            pts_indx = pts_indx([2 1]);
          end

          cortex = pos(ngroups < 5);
          if (cortex(1) < cortex(2))
            pts_indx = [fliplr(find(ngroups < 5)) pts_indx];
          else
            pts_indx = [find(ngroups < 5) pts_indx];
          end

          data(n, :) = data(n, pts_indx);

          disp([num2str(n) '/' num2str(nimg)]);
        end

        all_vals = [all_vals; [data ones(nimg, 1)*i]];
        all_noises = [all_noises; [noises ones(nimg, 1)*i]];
        disp([num2str(i) '/' num2str(nfiles)]);
      end

      noises_grp = mymean(all_noises, 1, all_noises(:,end));
      noise = noises_grp(all_vals(:,end),:);

      val = (all_vals(:,1) - all_vals(:,3)) ./ noise(:,2);
      m = mymean(val, 1, all_vals(:,end));
      val = (val - m(4)) / (m(5) - m(4));
      figure;boxplot(val, all_vals(:,end))
      title('SNR immuno')

      [H,p] = myttest(val, all_vals(:,end))


      targets = {'good_24.txt'; 'good_20.txt'; 'good_13.txt'; 'good_c27d91.txt'; 'good_ani2.txt'};
      data = load('data_expansion.mat');
      [junk, indxs] = ismember(targets, data.all_data(:,1));
      all_data = data.all_data(indxs,2);

      all_sizes = NaN(0,4);
      for i=1:length(targets)
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
      flows = [-advection_params advection_params(:,[end-1:-1:1])];
      figure;imagesc(flows)
      pos_indx = fliplr([size(advection_params,2):-25:1]);
      pos_indx = [pos_indx(1:end-1) pos_indx(end):25:size(advection_params,2)*2-1];
      pos = [1:2*size(advection_params,2)-1]-size(advection_params,2);

      t_indx = unique([fliplr([t0*10:-200:1]) t0*10:200:size(advection_params,1)]);

      set(gca, 'XTick', pos_indx, 'XTickLabel', pos(pos_indx))
      set(gca, 'YTick', t_indx, 'YTickLabel', t_indx - t0*10)
      colormap(redgreenmap)
      colorbar

      %figure;find_kymograph('1056-24-all.mat', 'config_fitting', 'fit_kymo', 'config_modeling', 'goehring', 'init_noise', 0, 'display', true, 'aligning_type', 'best');

      %figure;find_kymograph('1056-24-all.mat', 'config_fitting', 'fit_kymo', 'config_modeling', 'custom_flow', 'init_noise', 0, 'display', true, 'aligning_type', 'best');

      keyboard

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
      boxplot(all_times(:,1)-all_times(:,3), all_times(:,end), 'label', all_data(:,1));
      title('Polarity maintenance')
      [H,P] = myttest(all_times(:,1)-all_times(:,3), all_times(:,end))
      %boxplot((all_times(:,3)-all_times(:,2)) ./ (all_times(:,5)-all_times(:,4)), all_times(:,end), 'label', all_data(:,1));
      %title('Polarity establishment VS PC to PNM')
      %[H,P] = myttest((all_times(:,3)-all_times(:,2)) ./ (all_times(:,5)-all_times(:,4)), all_times(:,end))

      %subplot(2,3,6)
      %scatter((all_times(:,3)-all_times(:,2)), (all_times(:,5)-all_times(:,4)));
      %myregress([ones(size(all_times,1), 1) (all_times(:,3)-all_times(:,2))], (all_times(:,5)-all_times(:,4)));
      %title('Polarity establishment VS PC to PNM')

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
          %for t=1:length(thresh)
          %  tmp_indx(t) = find(fraction>=thresh(t), 1, 'first');
          %end
          tmp_indx = [find(fraction>=thresh(2), 1, 'first') length(fraction)];

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
        avgs = NaN(1,3);

      figure;hold on
      for i=1:3
        boxplot(all_data{i,2}(:,6), 'colors',colors(i, :), 'position', temperatures(i), 'width', 2);
        all_times = [all_times; [all_data{i,2}(:,6) ones(size(all_data{i,2}),1)*i]];

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

      title('Polarity maintenance')
      ylim([0 max(all_times(:,1)+200)]);
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
      opts = load_parameters(opts, 'extended_model.txt');

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
          timing = (length(path) - find(path>=thresh(2), 1, 'first'))*10;

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
      vals = group_ml_results('ScaledFlowFits/adr-kymo-*_evol.dat', {'parameter_set';'fit_flow';'fit_model';'fixed_parameter'}, {'type', '1056-temps-all'; 'aligning_type', 'fitting';'normalize_smooth', true; 'rescale_length_only', true; 'scale_each_egg', true; 'scale_flow', true});

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
      params = cell(size(score,1), 2);
      rel_params = params;
      orig_p = params;

      orig_v = vals;
      vals = extract_model_parameters(vals, true);

      for i=1:size(vals,1)
        ind = find(vals{i,1}{1}.fixed_parameter(5:end-1));
        if (isempty(ind))
          ind = 0;
        elseif (numel(ind) > 1)
          ind = 10*ind(1) + ind(2);
        end

        param_set(i,:) = [vals{i,1}{1}.parameter_set ind vals{i,1}{1}.fit_model vals{i,1}{1}.fixed_parameter(end)];
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

          params{i,1} = [vals{i,2}{indx,2}.params.rate vals{i,2}{indx,2}.params.offset vals{i,2}{indx,2}.params.energy vals{i,2}{indx,2}.params.viscosity vals{i,2}{indx,2}.params.flow vals{i,2}{indx,2}.params.sigma];
          rel_params{i,1} = [vals{i,2}{indx,2}.params.effective_value.viscosity; vals{i,2}{indx,2}.params.effective_value.rate; vals{i,2}{indx,2}.params.effective_value.flow];
          nparams(i) = vals{i,2}{indx,2}.params.nparams - sum(vals{i,1}{1}.fixed_parameter);
          params{i,2} = vals{i,1}{1}.fixed_parameter;
        else
          param_set(i,:) = NaN;
        end
      end

      [param_set, indx] = sortrows(param_set);
      goods = ~any(isnan(param_set),2);
      indx = indx(goods);
      score = score(indx);
      nparams = nparams(indx);
      params = params(indx,:);
      rel_params = rel_params(indx, 1);

      aic = 2*(nparams + score) + 2*nparams.*(nparams+1)./(npts-nparams-1);
      rel_L = exp((min(aic) - aic)/2)*100;

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
        vals = rel_params{i,1};
        vals = bsxfun(@rdivide, vals, vals(:,2));

        for j=2:size(vals,1)-1
          plot(temps, vals(j,:), '-s', 'Color', colors(j-1,:));
        end
        plot(temps, vals(end,:), '-v', 'Color', colors(end-1,:));
        plot(temps, vals(1,:), '-o', 'Color', colors(end,:));
        ylim([0 2.5]);
        xlim([10 26])

        title([num2str(param_set(i,:)) ':' num2str(aic(i)) ',' num2str(rel_L(i))]);
      end

      keyboard

    case 7.1
      vals = group_ml_results('LatestFits/adr-kymo-*_evol.dat', {'init_noise';'simulation_parameters';'init_pos'}, {'type', 'simulation'});
      colors = [37 82 115 150 zeros(1,10)] / 255;

      has_drawn_init_pos = false;

      is_log = [0 1 1 1 1 0];

      noise = NaN(size(vals,1),1);
      all_best = cell(size(noise));
      max_iter = 0;
      init_pos = [];
      all_init = NaN(0,4);

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
        if (~isempty(vals{i,1}{1}.init_pos))
          all_init(end+1,:) = vals{i,1}{1}.init_pos ./ vals{i,1}{1}.rescale_factor;
        end
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
          %pts(:,2:5) = abs(bsxfun(@rdivide, pts(:,2:5), vals{i,1}{1}.simulation_parameters ./ vals{i,1}{1}.rescale_factor));
          pts(:,6) = pts(:,6)*vals{i,1}{1}.offset_scaling;

          %best_vals(count,2:5) = abs(bsxfun(@rdivide, best_vals(count,2:5), vals{i,1}{1}.simulation_parameters ./ vals{i,1}{1}.rescale_factor));
          best_vals(count,6) = best_vals(count,6)*vals{i,1}{1}.offset_scaling;
          %init_pos(2:5) = init_pos(2:5) ./ vals{i,1}{1}.simulation_parameters;

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

      files = {'1056-ani2-290612_3_.mat', '1056-24-090311_11_.mat', '1056-c27d91-310712_1_.mat'};
      frames = [59, 64, 90];
      %files = {'1056-ani2-250612_3_.mat', '1056-24-250511_0_.mat', '1056-c27d91-040712_3_.mat'};
      %frames = [61, 63, 67];

      load(files{1});
      nimg = frames(1);

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
      title(files{1});

      figure;imagesc(domain)
      hold on;
      plot(path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      plot(-path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      set(gca, 'XTick', [-pos_indx([end:-1:2]) pos_indx]+width+1, 'XTickLabel', pos([-pos_indx([end:-1:2]) pos_indx]+width+1));
      set(gca, 'YTick', y_tick, 'YTickLabel', y_labels);
      title(files{1});
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
      title([files{1} ' normalized']);
      colormap(blueredmap)

      load(files{2});
      nimg = frames(2);

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
      title(files{2});

      figure;imagesc(domain)
      hold on;
      plot(path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      plot(-path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      set(gca, 'XTick', [-pos_indx([end:-1:2]) pos_indx]+width+1, 'XTickLabel', pos([-pos_indx([end:-1:2]) pos_indx]+width+1));
      set(gca, 'YTick', y_tick, 'YTickLabel', y_labels);
      title(files{2});
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
      title([files{2} ' normalized']);
      colormap(blueredmap)

      load(files{3});
      nimg = frames(3);

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
      title(files{3});

      figure;imagesc(domain)
      hold on;
      plot(path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      plot(-path+width+1, 1:length(path), 'Color', [83 83 83]/255);
      set(gca, 'XTick', [-pos_indx([end:-1:2]) pos_indx]+width+1, 'XTickLabel', pos([-pos_indx([end:-1:2]) pos_indx]+width+1));
      set(gca, 'YTick', y_tick, 'YTickLabel', y_labels);
      title(files{3});
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
      title([files{3} ' normalized']);
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

      %{
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
      %}

      avgs = NaN(3, 4);
      figure;

      subplot(2,3,6);
      myregress(all_params(:,2), all_params(:,1));
      title('Domain size vs membrane length, no intercept')
      xlabel('Membrane length (µm)');
      ylabel('Domain length (\mu m)')
      ylim([40 100]);
      xlim([80 200]);

      for i=1:3
        all_params = all_data{i,2};

        [means, stds] = mymean([all_params(:,2) all_params(:,1)]);

        subplot(2,3,i);hold on
        title(['N=' num2str(size(all_params,1))])
        xlabel('Membrane length (µm)');
        ylabel('Domain length (\mu m)')
        ylim([40 100]);
        xlim([80 200]);
        errorbarxy(means(1), means(2), stds(1), stds(2)); hold on
        myregress(all_params(:,2), all_params(:,1), colors(i, :));

        subplot(2,3,4);hold on
        ylim([40 100]);
        xlim([80 200]);
        errorbarxy(means(1), means(2), stds(1), stds(2)); hold on
        myregress(all_params(:,2), all_params(:,1), colors(i, :));
      end
      title('ani-2, wt, C27D9.1')
      xlabel('Egg length (µm)');
      ylabel('Fraction (a.u.)')
      ylim([40 100]);
      xlim([80 200]);

      keyboard

    case 4.2

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

      %colors = [254 217 166; ...
      %          251 180 174; ...
      %          222 203 228; ...
      %          255 127 0; ...
      %          228 26 28; ...
      %          152 78 163]/255;

      all_data = cell(3, 3);

      %data = load('data_expansion.mat');
      data = load('simul_expansion.mat');
      fits = load('data_fitting');
      %data.all_data = data.all_simul;

      targets = {'good_13.txt';'good_20.txt';'good_24.txt';'good_ani2.txt'; 'good_c27d91.txt'; 'averages'};
      [junk, indxs] = ismember(targets, fits.all_data(:,1));
      all_data = fits.all_data(indxs,:);
      sizes = cumsum(cellfun(@(x)(size(x,1)), all_data(:,2)));
      all_data = cat(1, all_data{1:end-1,2});
      all_data = [all_data(:,2:5) 2*all_data(:,end-3) all_data(:,1)./all_data(:,end)];
      outliers = pcout(all_data(:, 1:4), [0.33 25], [0.25 1]);

      targets = {'good_ani2.txt';'good_24.txt';'good_c27d91.txt'};

      nsims = size(data.all_simul, 1);

      x_abs = [80:10:200].';
      x_rel = [30:5:80].';
      x_rel = bsxfun(@power, x_rel, [0:1]);

      y_abs = NaN(length(x_abs), nsims);
      y_rel = NaN(size(x_rel, 1), nsims);

      all_pts = NaN(0, 3);

      warning off;

      figure;
      for p = 1:nsims
        %for f=1:size(data.all_simul{p}, 1)
          [junk, indxs] = ismember(targets, data.all_simul{p}(:,1));

          all_params = cat(1, data.all_simul{p}{indxs, 2});
          intercept = ones(size(all_params, 1), 1);

          b_abs = myregress(all_params(:,2), all_params(:,1));
          b_rel = myregress([intercept all_params(:,3)*2], all_params(:,1)./all_params(:,2));


          %subplot(2,3,1);hold on;
          %scatter(all_params(:,2), all_params(:,1));
          %subplot(2,3,2);hold on;
          %scatter(all_params(:,3), all_params(:,1)./all_params(:,2));
          subplot(2,3,1);hold on;
          if (outliers(p))
            y_tmp = x_abs*b_abs;
            plot(x_abs, y_tmp, '--', 'Color', colors(find(p<sizes, 1), :));
          else
            y_abs(:,p) = x_abs*b_abs;
            plot(x_abs, y_abs(:,p), 'Color', colors(find(p<sizes, 1), :));
          end
          subplot(2,3,2);hold on;
          if (outliers(p))
            y_tmp = sum(bsxfun(@times, x_rel, b_rel(:).'), 2);
            plot(x_rel(:,2), y_tmp, '--', 'Color', colors(find(p<sizes, 1), :));
          else
            y_rel(:,p) = sum(bsxfun(@times, x_rel, b_rel(:).'), 2);
            plot(x_rel(:,2), y_rel(:,p), 'Color', colors(find(p<sizes, 1), :));
          end

          all_pts = [all_pts; all_params(:,1:3)];

          %{
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
          %}

        %end
        disp([num2str(p) '/' num2str(nsims)])
      end
      %intercept = ones(size(all_pts, 1), 1);

      %subplot(2,3,5);hold on;
      %myregress(all_pts(:,2), all_pts(:,1));
      %ylim([40 100]);
      %xlim([80 200]);
      %subplot(2,3,6);hold on;
      %myregress([intercept all_pts(:,3)*2], all_pts(:,1)./all_pts(:,2));
      %ylim([0 1]);
      %xlim([30 80]);

      [m_abs, s_abs] = mymean(y_abs, 2);
      [m_rel, s_rel] = mymean(y_rel, 2);

      c_abs = x_abs\m_abs;
      c_rel = x_rel\m_rel;

      subplot(2,3,1);hold on;
      plot(x_abs, m_abs, 'k');
      errorbar(x_abs(2:2:end), m_abs(2:2:end), s_abs(2:2:end), 'k')
      title(c_abs)
      ylim([40 100]);
      xlim([80 200]);
      subplot(2,3,2);hold on;
      plot(x_rel(:,2), m_rel, 'k');
      errorbar(x_rel(2:2:end, 2), m_rel(2:2:end), s_rel(2:2:end), 'k')
      title(c_rel)
      ylim([0 1]);
      xlim([30 80]);

      warning on;

      keyboard

      %{
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
      %}

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

      tmp_all = cat(1,all_data{:,2});
      tmp_all = tmp_all(:,2:5);

      all_pts = NaN(0, 6);

      for i=1:length(targets)-1
        pts = all_data{i,2};
        pts = [pts(:,2:5) 2*pts(:,end-3) pts(:,1)./pts(:,end)];

        %pts = [pts(:,1).*pts(:,4) pts(:,2) pts(:,3).*pts(:,2) pts(:,4:end)];
        goods = ~outliers(1:size(pts,1), 1);
        outliers = outliers(size(pts,1)+1:end, 1);

        [means, stds] = mymean(pts(goods,:));

        if (i==1)
          hfig(1) = figure('units','normalized','outerposition',[0 0 1 1]);
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
          mtit(['Normalized score : ' num2str(outliers_fraction)]);
          hfig(2) = figure('units','normalized','outerposition',[0 0 1 1]);
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

              [rho, p] = corr(tmp_all(:,k), tmp_all(:,l));
              title(['r=' num2str(rho) ', p=' num2str(p)])
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
      vals = group_ml_results('FitsMarch/adr-kymo-*_evol.dat', {'combine_data'; 'fixed_parameter'}, {'type', '1056-all-all'; 'fitting_type', 'sample'; 'fit_relative', true});

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
      vals2 = vals2(cellfun(@(x)(length(x{1}.fixed_parameter) == 159), vals2(:,1)),:);

      if (strncmp(vals2{1,1}{1}.combine_data, 'hessian', 7))
        vals2 = vals2([2 1],:);
      end

      sensitivity_analysis(vals2{1,2}{1,2}.evolution, 'best_model');

      [rC, C, rH, H] = correlation_matrix(vals2{2,2}{1,2}.evolution);

      figure;imagesc(rC, [-1 1]);
      colormap(flipud(redgreencmap));
      colorbar

      keyboard

    case 8
      figs_msb(8.1);
      figs_msb(8.2);
    case 8.1
      colors = [254 217 166; ...
                251 180 174; ...
                222 203 228; ...
                255 127 0; ...
                228 26 28; ...
                152 78 163]/255;

      data = load('simul_optimized.mat');
      ref = load('data_expansion.mat');
      targets = {'good_ani2.txt';'good_24.txt';'good_c27d91.txt'};

      [junk, indxs] = ismember(targets, ref.all_data(:,1));
      all_ref = ref.all_data(indxs,2);
      sizes = cumsum(cellfun(@(x)(size(x,1)), all_ref));
      group_indx = ones(sizes(end), 1);
      group_indx(sizes(1)+1:sizes(2)) = 2;
      group_indx(sizes(2)+1:end) = 3;
      all_ref = cat(1, all_ref{:});
      all_ref = [2*all_ref(:,3) all_ref(:,1)./all_ref(:,2)];

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
        pval = chow_test([intercept all_params(:,6) all_params(:,1)./all_params(:,2)], [intercept all_ref]);
        title(pval)
        xlabel('Egg length (µm)');
        ylabel('Fraction (a.u.)')

        subplot(2,3,6);hold on;
        boxplot(all_params(:,1)./all_params(:,2), group_indx);
        ylim([0 1]);
        title(ranksum(all_ref(:,2), all_params(:,1)./all_params(:,2)));

        mtit('Simulation of membrane length fraction')


        avgs = NaN(3, 4);
        figure('Name', num2str(p));

        subplot(2,3,6);
        myregress(all_params(:,2), all_params(:,1), black, 'dx');
        title('Domain size vs membrane length, no intercept')
        xlabel('Membrane length (µm)');
        ylabel('Domain length (\mu m)')
        ylim([40 100]);
        xlim([80 200]);

        for i=1:3
          all_params = all_data{i,2};

          [means, stds] = mymean([all_params(:,2) all_params(:,1)]);

          subplot(2,3,i);hold on
          title(['N=' num2str(size(all_params,1))])
          xlabel('Membrane length (µm)');
          ylabel('Domain length (\mu m)')
          ylim([40 100]);
          xlim([80 200]);
          errorbarxy(means(1), means(2), stds(1), stds(2)); hold on
          myregress(all_params(:,2), all_params(:,1), colors(i, :), 'dx');

          subplot(2,3,4);hold on
          ylim([40 100]);
          xlim([80 200]);
          errorbarxy(means(1), means(2), stds(1), stds(2)); hold on
          myregress(all_params(:,2), all_params(:,1), colors(i, :), 'dx');
        end
        title('ani-2, wt, C27D9.1')
        xlabel('Egg length (µm)');
        ylabel('Fraction (a.u.)')
        ylim([40 100]);
        xlim([80 200]);


        [ratios, surfaces, volumes] = surface2volume(all_params(:,3:5).');
        heights = find_cap_size(all_params(:,3:5).', all_params(:,1).');
        [caps_surface, caps_volume] = spheroidal_cap(all_params(:,3:5).', heights);

        %{
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
        %}
      end

      keyboard

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
        indexes = cumsum(cellfun(@(x)(size(x,1)), data.all_simul{p}(:,2)));
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
        temperatures = data.temperatures(indexes(indxs));
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

    case 8.3

      colors = [251 180 174; ...
                204 235 197; ...
                179 205 227; ...
                222 203 228; ...
                254 217 166]/255;


      data = load('simul_optimized.mat');
      ref = load('data_expansion.mat');

      groups = [];

      for p=1:length(data.all_simul)
        all_ratios = NaN(0, 7);
        figure('Name', num2str(p));

        for f=1:size(data.all_simul{p}, 1)
          all_params = data.all_simul{p}{f,2};
          all_data = ref.all_data{f,2};

          all_values = [all_data(:,1) all_params(:,1) all_data(:,1)./all_data(:,2) all_params(:,1)./all_params(:,2)];
          %all_values = [all_params(:,1) all_data(:,1) all_params(:,1)./all_params(:,2) all_data(:,1)./all_data(:,2)];

          files = data.all_simul{p}{f,3};
          nfiles = size(files,1);
          for i=1:nfiles
            tmp_indx = NaN(2, length(thresh));
            fraction = files{i,2};
            ref_frac = ref.all_data{f,3}{i,2};
            %for t=1:length(thresh)
            %  tmp_indx(1, t) = find(fraction>=thresh(t), 1, 'first');
            %  tmp_indx(2, t) = find(ref_frac>=thresh(t), 1, 'first');
            %end
            tmp_indx(1, :) = [find(fraction>=thresh(2), 1, 'first') length(fraction)];
            tmp_indx(2, :) = [find(ref_frac>=thresh(2), 1, 'first') length(ref_frac)];

            all_values(i,5:6) = [(tmp_indx(2,2)-tmp_indx(2,1)) (tmp_indx(1,2)-tmp_indx(1,1))]*10;
          end

          intercept = ones(size(all_values,1), 1);
          [means, stds] = mymean(all_values);

          subplot(2,3,1);hold on
          errorbarxy(means(1), means(2), stds(1), stds(2)); hold on
          myregress([intercept all_values(:,1)], all_values(:,2), colors(f, :));
          text(0, 0, ref.all_data{f, 1}, 'Color', colors(f,:))
          ylim([30 80]);
          xlim([30 80]);

          subplot(2,3,2);hold on
          errorbarxy(means(3), means(4), stds(3), stds(4)); hold on
          myregress([intercept all_values(:,3)], all_values(:,4), colors(f, :));

          subplot(2,3,3);hold on
          errorbarxy(means(5), means(6), stds(5), stds(6)); hold on
          myregress([intercept all_values(:,5)], all_values(:,6), colors(f, :));
          ylim([0 1500]);
          xlim([0 1500]);

          all_ratios = [all_ratios; [all_values intercept*f]];
        end

        intercept = ones(size(all_ratios,1), 1);
        [means, stds] = mymean(all_ratios);

        subplot(2,3,4);hold on
        errorbarxy(means(1), means(2), stds(1), stds(2)); hold on
        myregress([intercept all_ratios(:,1)], all_ratios(:,2));
        [rho, p] = corr(all_ratios(:,1), all_ratios(:,2));
        title(['r=' num2str(rho) ', p=' num2str(p)])
        ylim([30 80]);
        xlim([30 80]);

        subplot(2,3,5);hold on
        errorbarxy(means(3), means(4), stds(3), stds(4)); hold on
        myregress([intercept all_ratios(:,3)], all_ratios(:,4));
        [rho, p] = corr(all_ratios(:,3), all_ratios(:,4));
        title(['r=' num2str(rho) ', p=' num2str(p)])

        subplot(2,3,6);hold on
        errorbarxy(means(5), means(6), stds(5), stds(6)); hold on
        myregress([intercept all_ratios(:,5)], all_ratios(:,6), [0 0 1]);
        [rho, p] = corr(all_ratios(:,5), all_ratios(:,6));
        title(['r=' num2str(rho) ', p=' num2str(p)])
        ylim([0 1500]);
        xlim([0 1500]);
      end

      keyboard

    case 9

      load('data_expansion');
      tmp_data = cat(1, all_data{:,3});
      all_names = tmp_data(:, 3);

      tmp_opts = get_struct('modeling');
      opts_goehring = load_parameters(tmp_opts, 'goehring.txt');
      opts_flow = load_parameters(opts_goehring, 'custom_flow.txt');
      opts_size = load_parameters(opts_goehring, 'extended_model.txt');
      opts_full = load_parameters(opts_goehring, 'final_model.txt');

      flow_goehring = size(opts_goehring.advection_params, 2);

      orig_flow = [opts_flow.diffusion_params; opts_flow.reaction_params];
      orig_flow([4 12]) = orig_flow([4 12]) ./ orig_flow([13 5]);
      orig_flow = orig_flow(:).';

      orig_size = [opts_size.diffusion_params; opts_size.reaction_params];
      orig_size([4 12]) = orig_size([4 12]) ./ orig_size([13 5]);
      orig_size = orig_size(:).';

      orig_full = [opts_full.diffusion_params; opts_full.reaction_params];
      orig_full([4 12]) = orig_full([4 12]) ./ orig_full([13 5]);
      orig_full = orig_full(:).';

      all_scores = Inf(length(all_names), 5);
      all_offsets = all_scores;

      vals = [group_ml_results('LatestFits/adr-kymo-*_evol.dat', {'parameter_set', 2; 'scale_flow', false; 'fit_flow', false}, {'type';'simulation_parameters';'fixed_parameter'; 'flow_size'}); ...
              group_ml_results('LatestFits/adr-kymo-*_evol.dat', {'parameter_set', 0; 'scale_flow', false; 'fit_flow', false}, {'type';'simulation_parameters';'fixed_parameter';'flow_size'}); ...
              group_ml_results('FitsMarch/adr-kymo-*_evol.dat', {'parameter_set', 0; 'scale_flow', true}, {'type';'simulation_parameters';'fixed_parameter';'fit_flow'})];
      vals = extract_model_parameters(vals, true);

      nfits = size(vals,1);

      for i=1:nfits
        curr_val = vals{i,1}{1};

        sub_indx = 0;

        indx = find(ismember(all_names, [curr_val.type '.mat']), 1);
        if (isempty(indx))
          continue;
        end

        if (curr_val.parameter_set==2)
          sub_indx = 5;
        else
          if (flow_goehring == curr_val.flow_size(2))
            sub_indx = 1;
          elseif (numel(curr_val.simulation_parameters)>4 && sum((curr_val.simulation_parameters(4:5) - orig_flow(4:5)).^2)<1e-6)
              sub_indx = 2;
          elseif (numel(curr_val.simulation_parameters)>4 && sum((curr_val.simulation_parameters(4:5) - orig_size(4:5)).^2)<1e-6)
              sub_indx = 3;
          elseif (numel(curr_val.simulation_parameters)>4 && sum((curr_val.simulation_parameters(4:5) - orig_full(4:5)).^2)<1e-6)
              sub_indx = 4;
          end
        end

        if (sub_indx ~= 0 && all_scores(indx, sub_indx) > vals{i,2}{1,2}.score)
          all_scores(indx, sub_indx) = vals{i,2}{1,2}.score;
          all_offsets(indx, sub_indx) = vals{i,2}{1,2}.params.offset;
        end
      end

      figure;bar(sum(all_scores))

      keyboard

      %{
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
      %}

    case 10
      vals = group_ml_results('LatestFits/adr-kymo-*_evol.dat', {'parameter_set';'fitting_type'; 'scale_flow'; 'extrapol_z'}, {'type', '1056-med-scale'});
      vals = [vals; group_ml_results('LatestFits/adr-kymo-*_evol.dat', {'parameter_set';'fitting_type'; 'scale_flow'; 'extrapol_z'}, {'type', '1056-med-all'})];
      values = extract_model_parameters(vals, true);

      all_scales = cell(0,2);
      count = 1;

      for i=1:size(values, 1)
        for j=1:size(values{i,2}, 1)
          params = [values{i,2}{j,2}.params.rate values{i,2}{j,2}.params.offset values{i,2}{j,2}.params.energy values{i,2}{j,2}.params.viscosity values{i,2}{j,2}.params.flow values{i,2}{j,2}.params.sigma values{i,2}{j,2}.params.flow_scaling];

          [score, results] = check_parameters([values{i,1}{1}.type '.mat'], values{i,1}{1}, 'config_modeling', 'custom_flow', 'init_pos', params, 'start_with_best', false);
          tmp_results = NaN(length(results), 2);

          for g = 1:length(results)
            domain = results{g, 1};

            domain = domain((end/2)+1:end, :).';
            domain = [domain domain(:,end-1:-1:1)];

            [profile, center, max_width, cell_width, path] = get_profile(domain, nframes);

            nprofile = length(profile)-1;
            dwidth = mymean(results{g,2}(end-nframes+1:end));

            tmp_results(g, :) = [score, 2*dwidth];
          end

          all_scales{count, 1} =  tmp_results;
          all_scales{count, 2} =  results;

          count = count + 1;
        end
      end

      keyboard

    case 20
      vals = group_ml_results('FitsRevision/adr-kymo-*_evol.dat', {'parameter_set',2});
      values = extract_model_parameters(vals, true);

      keyboard

    case 21
      x = [0:2:14];
      y = [-0.2:0.1:0.3];

      xp = [169:69:660];
      yp = [412 345 278 211 144 77];

      data = [194, 346; 245, 322; 298, 291; 348, 249; 401, 213; 452 200; 504, 187; 555, 179; 607, 178];
      stds = [194, 321; 245 346; 298, 315; 348, 224; 401, 193; 452, 180; 504, 168; 555, 156; 607, 153];

      vals = [interp1(xp, x, data(:,1).'); interp1(yp, y, data(:,2).')].';
      stds = [interp1(xp, x, stds(:,1).'); interp1(yp, y, stds(:,2).')].';
      stds = abs(vals(:,2) - stds(:,2));

      load('niwayama.mat')

      %nparticles = 2000;
      for nparticles = [100:200:1000]
        xr = rand([nparticles, 1])*range(x)+x(1);
        yr = randn(nparticles, 1).*feval(fittedstds, xr) + feval(fittedmodel, xr);

        edges = [0:1.5:14];
        centers = edges(1:end-1) + diff(edges)/2;

        edges(1) = -1e-3;
        edges(end) = 14;

        [n, bin] = histc(xr, edges);

        newvals = vals;
        for i=1:size(vals,1)
          [newvals(i,1), newvals(i,2)] = mymean(yr(bin==i));
        end

        pos = [x(1):0.1:x(end)];
        spline = feval(fittedmodel, pos);
        splinestds = feval(fittedstds, pos);

        weight = @(x)(exp(-(x.^2)/2));
        weights = weight(xr);
        weights = weights / sum(weights);

        avg = sum(yr .* weights);
        poss = [spline(1) vals(1,2) avg];

        figure;scatter(xr, yr);hold on
        errorbar(vals(:,1), vals(:,2), stds, 'r')
        errorbar(centers, newvals(:,1), newvals(:,2), 'g')
        plot(pos, spline+splinestds, 'k')
        plot(pos, spline-splinestds, 'k')
        plot(pos, spline, 'k')
        plot(pos, weight(pos)*max(yr), 'm')
        scatter(pos([1 1]), [spline(1) avg], 'r')
        title([num2str(nparticles) ':' num2str(poss) ' ;' num2str(abs(poss - poss(1)))])
      end

      keyboard

    case 22

      %{
      load('749-031012_5_');

      all_movies = cell(4,2);
      for resc = 1:4
        for rem = 0:1
          new_opts = opts;
          new_opts.spot_tracking.projection_args = [Inf resc rem];
          all_movies{resc, rem+1} = piv_flows(mymovie, new_opts);

          disp([num2str(rem) ',' num2str(resc)]);
        end
      end

      save([mymovie.experiment 'pivs'], 'all_movies', 'opts');

      all_flows = cell(size(all_movies));
      for i=1:size(all_movies,1)
        for j=1:size(all_movies,2)
          all_flows{i,j}=display_flow(all_movies{i,j}, opts, true);
        end
      end
      %}

      load('749-031012_5_pivs')

      resc=2;
      rem=0;
      %clim = [-0.5 0.5];
      clim = [-0.16 0.16];

      mymovie = all_movies{resc, rem+1};
      [orig_flow, pos] = display_flow(mymovie, opts, false);
      piv_flow = display_flow(mymovie, opts, true);

      marks = [-50:25:50];
      real_pos = interp1(pos.', 1:length(pos), marks);

      figure;subplot(1,2,1);
      colormap(redgreenmap)
      imagesc(orig_flow.', clim);
      set(gca, 'XTick', real_pos, 'XTickLabel', marks)
      subplot(1,2,2);
      imagesc(piv_flow.', clim);
      set(gca, 'XTick', real_pos, 'XTickLabel', marks)
      colorbar

      keyboard

      t0 = 35; % corr_offset + yindx in the next find_kymograph

      advection_params = getfield(load('cyto_flow.mat', 'flow3'), 'flow3');
      flows = [-advection_params advection_params(:,[end-1:-1:1])];
      figure;imagesc(flows)
      pos_indx = fliplr([size(advection_params,2):-25:1]);
      pos_indx = [pos_indx(1:end-1) pos_indx(end):25:size(advection_params,2)*2-1];
      pos = [1:2*size(advection_params,2)-1]-size(advection_params,2);

      t_indx = unique([fliplr([t0*10:-200:1]) t0*10:200:size(advection_params,1)]);

      set(gca, 'XTick', pos_indx, 'XTickLabel', pos(pos_indx))
      set(gca, 'YTick', t_indx, 'YTickLabel', t_indx - t0*10)
      colormap(redgreenmap)
      colorbar



      keyboard

    case 23

      clim = [-0.16 0.16];
      resc=2;
      rem=0;
      files = dir('749-*_.mat');
      all_flows = cell(length(files), 2);
      align = load('aligned_flow_1_1127907.619.mat');

      t_ref = 457;

      for i=1:length(files)
        load(files(i).name);

        if (~isfield(mymovie.data, 'piv'))
          new_opts = opts;
          new_opts.spot_tracking.projection_args = [Inf resc rem];
          mymovie = piv_flows(mymovie, new_opts);
          save(mymovie.experiment, 'mymovie', 'opts');
        end

        [orig_flow, pos] = display_flow(mymovie, opts, false);
        piv_flow = display_flow(mymovie, opts, true);
        piv_flow(:,1:12) = NaN;
        switch i
          case 4
            piv_flow(:,[13:25]) = NaN;
          case 6
            piv_flow(:,[65:75]) = NaN;
        end

        all_flows{i, 1} = orig_flow.';
        all_flows{i, 2} = piv_flow.';

        if (i==1)
          marks = [-50:25:50];
          real_pos = interp1(pos.', 1:length(pos), marks);
          t0 = t_ref - (max(align.signals_full{1,2}) - align.signals_full{1,2}(i));
          t_indx = unique([fliplr([t0:-200:1]) t0:200:size(orig_flow,2)]);

          figure;subplot(1,2,1);
          colormap(redgreenmap)
          imagesc(orig_flow.', clim);
          set(gca, 'XTick', real_pos, 'XTickLabel', marks)
          set(gca, 'YTick', t_indx, 'YTickLabel', t_indx - t0)
          subplot(1,2,2);
          imagesc(piv_flow.', clim);
          set(gca, 'XTick', real_pos, 'XTickLabel', marks)
          set(gca, 'YTick', t_indx, 'YTickLabel', t_indx - t0)
          colorbar
        end

        disp([num2str(i) '/' num2str(length(files))]);
      end

      avg_flow = average_flow(all_flows(:,1), align.signals_full{1,2}, 2);
      avg_piv = average_flow(all_flows(:,2), align.signals_full{1,2}, 2);

%      figure;subplot(1,2,1);
%      colormap(redgreenmap)
%      imagesc(avg_flow(:,:,3), clim);
      %set(gca, 'XTick', real_pos, 'XTickLabel', marks)
%      subplot(1,2,2);
%      imagesc(avg_piv(:,:,3), clim);
      %set(gca, 'XTick', real_pos, 'XTickLabel', marks)
      %colorbar

      t0 = 35; % corr_offset + yindx in the next find_kymograph

      advection_params = avg_flow(:,:,3);
      [m,s] = mymean(advection_params(advection_params > mean(advection_params(:)+2*std(advection_params(:)))));
      flows = [-advection_params advection_params(:,[end-1:-1:1])];
      figure;imagesc(flows)
      pos_indx = fliplr([size(advection_params,2):-25:1]);
      pos_indx = [pos_indx(1:end-1) pos_indx(end):25:size(advection_params,2)*2-1];
      pos = [1:2*size(advection_params,2)-1]-size(advection_params,2);

      t_indx = unique([fliplr([t0*10:-200:1]) t0*10:200:size(advection_params,1)]);

      set(gca, 'XTick', pos_indx, 'XTickLabel', pos(pos_indx))
      set(gca, 'YTick', t_indx, 'YTickLabel', t_indx - t0*10)
      title([num2str(m) '+-' num2str(s)])
      colormap(redgreenmap)
      colorbar

      advection_params = avg_piv(:,:,3);
      [m,s] = mymean(advection_params(advection_params > mean(advection_params(:)+2*std(advection_params(:)))));
      flows = [-advection_params advection_params(:,[end-1:-1:1])];
      figure;imagesc(flows)
      pos_indx = fliplr([size(advection_params,2):-25:1]);
      pos_indx = [pos_indx(1:end-1) pos_indx(end):25:size(advection_params,2)*2-1];
      pos = [1:2*size(advection_params,2)-1]-size(advection_params,2);

      t_indx = unique([fliplr([t0*10:-200:1]) t0*10:200:size(advection_params,1)]);

      set(gca, 'XTick', pos_indx, 'XTickLabel', pos(pos_indx))
      set(gca, 'YTick', t_indx, 'YTickLabel', t_indx - t0*10)
      title([num2str(m) '+-' num2str(s)])
      colormap(redgreenmap)
      colorbar

      all_vals = all_flows(:,1);
      all_vals = cat(1, all_vals{:});
      all_vals = abs(all_vals(:));

      all_pivs = all_flows(:,2);
      all_pivs = cat(1, all_pivs{:});
      all_pivs = abs(all_pivs(:));

      figure;
      distributionPlot([all_vals;all_pivs], 'groups', [zeros(size(all_vals));ones(size(all_pivs))],'histOpt', 1.1,'xNames', {'Particles', 'PIV'});
      ylim([0 0.25])

      keyboard

    case 24

      filesw = [dir('749-0*_.mat');dir('749-2*_.mat')];
      filesa = dir('749-a*_.mat');
      filesa = filesa(2:end);
      %filesa = filesa([2 4:end]);
      filesc = dir('749-c*_.mat');
      %filesc = filesc(2:end);
      %filesc = filesc(2:end-1);
      files = [filesw; filesa; filesc];
      %files = [filesa; filesc];
      nfiles = length(files);
      nwt = length(filesw);
      %nwt = 0;
      nani = length(filesa);
      nc27 = length(filesc);
      all_flows = cell(nfiles, 1);
      all_sizes = NaN(3,nfiles);
      all_stats = NaN(nfiles,10);

      convert = get_struct('z-correlation');
      names = cell(nfiles, 1);

      for i=1:nfiles
        load(files(i).name);
        names{i} = files(i).name;
        orig_ps = opts.pixel_size;
        ps_ratio = 8/16.51;
        if (i>nwt)
          opts.ccd_pixel_size = 16.51;
          opts = set_pixel_size(opts);
        else
          opts.ccd_pixel_size = 8;
          opts = set_pixel_size(opts);
        end
        all_flows{i} = display_flow(mymovie, opts);
        all_flows{i} = (all_flows{i}(1:end/2,:) - all_flows{i}(end/2+1:end,:))/2;
        if (i>nwt)
          all_flows{i} = abs(all_flows{i}(:))*ps_ratio;
          %all_flows{i} = abs(all_flows{i}(:))*opts.pixel_size;
          %all_flows{i} = abs(all_flows{i}(:))*(opts.pixel_size/orig_ps);
        else
          all_flows{i} = abs(all_flows{i}(:));
          %all_flows{i} = abs(all_flows{i}(:))*(opts.pixel_size/orig_ps);
        end
        cortex = mymovie.data.cortex(end).carth;
        conv = convhull(cortex(:,1),cortex(:,2));
        [junk,all_sizes(1:2,i),junk] = fit_ellipse(cortex(conv,:));
        all_sizes(:,i) = all_sizes(:,i)*opts.pixel_size;

        all_stats(i,1:8) = [prctile(all_flows{i},[50:10:100]) nanmean(all_flows{i}) nanstd(all_flows{i})];
        [all_stats(i,9), all_stats(i,10)] = mymean(all_flows{i}(all_flows{i}>=all_stats(i,7) + 2*all_stats(i,8)));

        disp([num2str(i) '/' num2str(nfiles)]);
      end

      %all_sizes = [[29;19.5;15], all_sizes];

      all_sizes(3, :) = convert.bkg + convert.long_axis*all_sizes(1, :) + convert.short_axis*all_sizes(2, :);

      all_sizes2 = bsxfun(@times, all_sizes, [1.034; 1.016; 1]);
      all_sizes2(3, :) = convert.bkg + convert.long_axis*all_sizes2(1, :) + convert.short_axis*all_sizes2(2, :);

      all_ratios = [surface2volume(all_sizes);surface2volume(all_sizes2)];

      ref_flow = load('cyto_flow.mat');
      ref_flow = ref_flow.flow3;
      ref_flow = abs(ref_flow(:));
      ref_stats=NaN(1,10);
      ref_stats(1,1:8) = [prctile(ref_flow,[50:10:100]) nanmean(ref_flow) nanstd(ref_flow)];
      [ref_stats(1,9), ref_stats(1,10)] = mymean(ref_flow(ref_flow>=ref_stats(1,7) + 2*ref_stats(1,8)));

      %ref = surface2volume([58;38;30]/2);
      %all_ratios(2,1) = all_ratios(1,1);
      %ref = all_ratios(1,1);
      ref = 1;

      %all_stats = [ref_stats; all_stats];

      %nfiles = nfiles + 1;
      %nwt = nwt + 1;
      %names = ['Ref'; names];

      g = 1.46;
      %preds = (1 - g*(ref./all_ratios - 1));
      preds = (1 + g*(ref./all_ratios - 1));
      intercept = ones(nfiles, 1);

      %for i=1:length(ref_stats)
      for i=9
        figure;

        subplot(2,2,1);hold on;
        myregress([intercept 2*all_sizes2(1,:).'], all_stats(:,i), 'b');
        errorbar(2*all_sizes2(1,:), all_stats(:,i), all_stats(:,i+1),'k','LineStyle', 'none');
        scatter(2*all_sizes2(1,1:nwt), all_stats(1:nwt,i),'g');
        scatter(2*all_sizes2(1,nwt+1:nwt+nani), all_stats(nwt+1:nwt+nani,i),'r');
        scatter(2*all_sizes2(1,nwt+nani+1:end), all_stats(nwt+nani+1:end,i),'m');
        %text(all_sizes2(1,:), all_stats(:,i),names);
        xlim([30 80])
        ylim([0 0.2])

        subplot(2,2,2);hold on;
        %scatter(all_sizes(1,1:6), all_stats(1:6,i),'b');
        %myregress([intercept (ref./all_ratios(2,:).') - 1], all_stats(:,i) - all_stats(1,i), 'b');
        myregress([(ref./all_ratios(2,:).') - 1], all_stats(:,i) - all_stats(1,i), 'b');
        errorbar((ref./all_ratios(2,:))-1, all_stats(:,i) - all_stats(1,i), all_stats(:,i+1),'k','LineStyle', 'none');
        scatter((ref./all_ratios(2,1:nwt))-1, all_stats(1:nwt,i) - all_stats(1,i),'g');
        scatter((ref./all_ratios(2,nwt+1:nwt+nani))-1, all_stats(nwt+1:nwt+nani,i) - all_stats(1,i),'r');
        scatter((ref./all_ratios(2,nwt+nani+1:end))-1, all_stats(nwt+nani+1:end,i) - all_stats(1,i),'m');
        %scatter((ref./all_ratios(2,:))-1, preds(1,:)*ref_stats(i),'k');
        text((ref./all_ratios(2,:))-1, all_stats(:,i) - all_stats(1,i),names);

        subplot(2,2,3);hold on;
        %scatter(all_sizes(1,1:6), all_stats(1:6,i),'b');
        myregress([intercept all_ratios(2,:).'], all_stats(:,i), 'b');
        errorbar(all_ratios(2,:), all_stats(:,i), all_stats(:,i+1),'k','LineStyle', 'none');
        scatter(all_ratios(2,1:nwt), all_stats(1:nwt,i),'g');
        scatter(all_ratios(2,nwt+1:nwt+nani), all_stats(nwt+1:nwt+nani,i),'r');
        scatter(all_ratios(2,nwt+nani+1:end), all_stats(nwt+nani+1:end,i),'m');
        %scatter((ref./all_ratios(2,:))-1, preds(1,:)*ref_stats(i),'k');
        %text(all_ratios(2,:), all_stats(:,i),names);
        xlim([0.16 0.18])
        ylim([0 0.2])

        subplot(2,2,4);hold on;
        %scatter(all_sizes(1,1:6), all_stats(1:6,i),'b');
        myregress([intercept 1./all_ratios(2,:).'], all_stats(:,i), 'b');
        errorbar(1./all_ratios(2,:), all_stats(:,i), all_stats(:,i+1),'k','LineStyle', 'none');
        scatter(1./all_ratios(2,1:nwt), all_stats(1:nwt,i),'g');
        scatter(1./all_ratios(2,nwt+1:nwt+nani), all_stats(nwt+1:nwt+nani,i),'r');
        scatter(1./all_ratios(2,nwt+nani+1:end), all_stats(nwt+nani+1:end,i),'m');
        %scatter((ref./all_ratios(2,:))-1, preds(1,:)*ref_stats(i),'k');
        text(1./all_ratios(2,:), all_stats(:,i),names);


        %subplot(2,2,3);hold on;
        %scatter(all_stats(:,i), preds(1,:)*ref_stats(i),'b');
        %myregress([intercept all_stats(:,i)], (preds(1,:).')*ref_stats(i), 'b');
        %subplot(2,2,4);hold on;
        %scatter(all_stats(:,i), preds(2,:)*ref_stats(i),'b');
        %myregress([intercept all_stats(:,i)], (preds(2,:).')*ref_stats(i), 'b');

      end

      keyboard

    case 25

      files = dir('1056-z-size_*_.mat');
      nfiles = length(files);
      sizes = NaN(3,nfiles);
      intercept = ones(nfiles, 1);

      for i=1:nfiles
        load(files(i).name);
        sizes(:,i) = mymovie.metadata.axes_length_3d;
      end

      convert = get_struct('z-correlation');

      figure;
      subplot(2,2,1);hold on
      myregress([intercept, sizes(1,:).'],sizes(3,:).');

      subplot(2,2,2);hold on
      myregress([intercept, sizes(2,:).'],sizes(3,:).');

      subplot(2,2,3);hold on
      myregress([intercept, sizes(1:2,:).'],sizes(3,:).');

      keyboard

    case 11
      %vals = group_ml_results('ScaledFlowFits/adr-kymo-*_evol.dat', {'type', '1056-24-all'}, {'simulation_parameters';'flow_size';'init_pos'});

      score = NaN(0,1);
      results = cell(0,2);

      [score(end+1), results(end+1,:)] = check_parameters('1056-24-all.mat', 'config_fitting', 'fit_kymo', 'config_modeling', 'goehring', 'aligning_type', 'best', 'start_with_best', false);
      [score(end+1), results(end+1,:)] = check_parameters('1056-24-all.mat', 'config_fitting', 'fit_kymo', 'config_modeling', 'custom_flow', 'aligning_type', 'best', 'start_with_best', false);
      [score(end+1), results(end+1,:)] = check_parameters('1056-24-all.mat', 'config_fitting', 'fit_kymo', 'config_modeling', 'custom_flow', 'aligning_type', 'best', 'start_with_best', false, 'init_pos', [0.00769 2.197 0.0314 2.202]);
      [score(end+1), results(end+1,:)] = check_parameters('1056-24-all.mat', 'config_fitting', 'fit_kymo', 'config_modeling', 'custom_flow', 'aligning_type', 'best', 'start_with_best', false, 'init_pos', [0.0116 2.1571 0.0658 2.1871]);
      [score(end+1), results(end+1,:)] = check_parameters('1056-24-all.mat', 'config_fitting', 'fit_kymo', 'config_modeling', 'extended_model', 'aligning_type', 'best', 'start_with_best', false);
      [score(end+1), results(end+1,:)] = check_parameters('1056-24-all.mat', 'config_fitting', 'fit_flows', 'config_modeling', 'custom_flow', 'aligning_type', 'best', 'start_with_best', false, 'init_pos', [0.002684 2.1571 0.0151 2.1871 0.3045 0 0.8478 0.8772 1.4575], 'parameter_set', 15);

      load('1056-24-all.mat');
      nlayers = size(ground_truth, 3);
      tmp_plot_truth = mymean(ground_truth,3).';

      boundary = (size(ground_truth,2)-1)/2;
      ground_truth = permute([ground_truth(:, 1:boundary+1, :), ground_truth(:,end:-1:boundary+1, :)], [2 1 3]);
      pos = pos(boundary+1:end);

      opts_expansion = load_parameters(get_struct('ASSET'), 'domain_expansion.txt');

      [f, frac_width, full_width] = domain_expansion(mymean(ground_truth(1:end/2, :, :), 3).', mymean(ground_truth((end/2)+1:end, :, :), 3).', size(ground_truth,1)/2, size(ground_truth,2), opts_expansion);
      tmp_gfrac = f*frac_width;
      tmp_gfrac(isnan(tmp_gfrac)) = 0;
      gfraction = tmp_gfrac;

      max_err = NaN(length(score), 1);
      residuals = cell(length(score), 1);

      for i=1:length(score)
        res = results{i,1};
        fraction = results{i,2};
        tmp_err = (ground_truth - repmat(res, [1 1 nlayers])).^2;

        tmp_plot_avg = res.';
        tmp_plot_err = mymean(tmp_err,3).';

        residuals{i} = tmp_plot_err;
        max_err(i) = max(tmp_plot_err(:));

        tmp_plot_frac = fraction;
        tmp_plot_frac(1:find(tmp_plot_frac, 1, 'first')-1) = NaN;
        tmp_plot_gfrac = gfraction;
        tmp_plot_gfrac(1:find(tmp_plot_gfrac, 1, 'first')-1) = NaN;

        yindx = find(tmp_plot_gfrac>0, 1, 'first');
        y_tick = unique([fliplr([yindx:-20:1]) yindx:20:size(tmp_plot_avg,1)]);
        y_labels = (y_tick - yindx)*10;

        xfrac = [(size(tmp_plot_err, 2)/2)-2*tmp_plot_frac(end:-1:1); (size(tmp_plot_err, 2)/2)+2*tmp_plot_frac];
        xgfrac = [(size(tmp_plot_err, 2)/2)-2*tmp_plot_gfrac(end:-1:1); (size(tmp_plot_err, 2)/2)+2*tmp_plot_gfrac];
        yfrac = [size(tmp_plot_err, 1):-1:1 1:size(tmp_plot_err,1)].';

        tmp_pos_tick = [0:50:(size(tmp_plot_err,2)/2 - 1)];
        tmp_pos_label = pos(tmp_pos_tick + 1);
        tmp_pos_tick = [-tmp_pos_tick(end:-1:2) tmp_pos_tick] + (size(tmp_plot_err,2)/2);
        tmp_pos_label = [-tmp_pos_label(end:-1:2) tmp_pos_label];

        figure;
        hax = subplot(1,3,1, 'NextPlot', 'replace');
        %hold off;
        imagesc([tmp_plot_avg(:,1:end/2) tmp_plot_avg(:,end:-1:(end/2)+1)], 'Parent', hax);
        set(hax, 'XTick', tmp_pos_tick, 'XTickLabel', tmp_pos_label, ...
                 'YTick', y_tick, 'YTickLabel', y_labels, ...
                 'NextPlot', 'add')
        colormap(hax, blueredmap)
        plot(hax, xfrac, yfrac, 'Color', [83 83 83]/255);

        hax = subplot(1,3,2, 'NextPlot', 'replace');
        %hold off;
        %imagesc(mymean(tmp_err, 3));
        imagesc([tmp_plot_err(:,1:end/2) tmp_plot_err(:,end:-1:(end/2)+1)], 'Parent', hax);
        set(hax, 'XTick', tmp_pos_tick, 'XTickLabel', tmp_pos_label, ...
                 'YTick', y_tick, 'YTickLabel', y_labels, ...
                 'NextPlot', 'add')
        colormap(hax, redblackmap)

        title(score(i));
        %hold on;
        %plot((size(tmp_err, 2)/2)-2*gfraction{g}, 1:size(tmp_err,1), 'k');
        %plot((size(tmp_err, 2)/2)+2*gfraction{g}, 1:size(tmp_err,1), 'k');
        %plot((size(tmp_err, 2)/2)-2*fraction, 1:size(tmp_err,1), 'w');
        %plot((size(tmp_err, 2)/2)+2*fraction, 1:size(tmp_err,1), 'w');
        plot(hax, xfrac, yfrac, 'Color', [83 83 83]/255);
        plot(hax, xgfrac, yfrac, 'Color', [83 83 83]/255);

        hax = subplot(1,3,3, 'NextPlot', 'replace');
        %hold off;
        imagesc(tmp_plot_truth.', 'Parent', hax);
        %set(gca, 'XTick', tmp_pos_tick, 'XTickLabel', tmp_pos_label);
        %set(gca, 'YTick', y_tick, 'YTickLabel', y_labels);
        set(hax, 'XTick', tmp_pos_tick, 'XTickLabel', tmp_pos_label, ...
                 'YTick', y_tick, 'YTickLabel', y_labels, ...
                 'NextPlot', 'add')
        colormap(hax, blueredmap)
        %hold on
        plot(hax, xgfrac, yfrac, 'Color', [83 83 83]/255);
      end

      for i=1:length(max_err)
        figure;
        colormap(redblackmap);
        for j=1:length(residuals);
          subplot(2,3,j);
          imagesc([residuals{j}(:,1:end/2) residuals{j}(:,end:-1:(end/2)+1)], [0 max_err(i)]);
        end
      end

      keyboard

    case 12

      best_align = load('1056-24-all.mat');
      files = textread('good_24.txt', '%s');
      segment = load('data_expansion');

      nfiles = length(files);
      prct_thresh = 5;

      if (~exist('all_kymos.mat'))
        all_files = cell(nfiles, 2);
        timings = NaN(nfiles, 3);
        scaling_factors = NaN(nfiles, 2);
        for i=1:nfiles
          data = load(files{i});
          [timings(i,:), name] = get_manual_timing(data.mymovie, data.opts);
          [f, frac_width, full_width, domain, pos] = domain_expansion(data.mymovie, data.opts);
          scaling_factors(i,1) = (size(domain, 2) - 1) / 2;
          scaling_factors(i,2) = frac_width/opts_expansion.quantification.resolution;

          all_files{i,1} = domain;
          all_files{i,2} = f;

          img = all_files{i,1};
          path = all_files{i,2}*scaling_factors(i, 2);
          h = scaling_factors(i,1);
          w = size(img, 1);
          pos_mat = repmat([1:h], w, 1);
          mask = bsxfun(@le, pos_mat, path);
          mask = [mask(:,end:-1:2) mask];

          min_val = prctile(img(~mask), prct_thresh);
          max_val = prctile(img(mask), 100-prct_thresh);
          if (min_val > max_val)
            all_files{i,1} = (img - max_val) / (min_val - max_val);
          else
            all_files{i,1} = (img - min_val) / (max_val - min_val);
          end

          disp([num2str(i) '/' num2str(nfiles)])
        end
        save('all_kymos', 'all_files', 'timings', 'scaling_factors')
      else
        load('all_kymos.mat')
      end

      dt = timings(:,2) - timings(:,1);
      dt2 = timings(:,3) - timings(:,1);
      ratio = nanmedian(dt ./ dt2);
      timings(isnan(dt),2) = timings(isnan(dt),1)+round(dt2(isnan(dt))*ratio);

      thresh = [0 5 10 25 50 75 77.5 90 100]/100;
      img_thresh = [0:0.1:0.9];

      msr = NaN(length(thresh)+3,length(img_thresh));
      for j=1:length(img_thresh)
        for i=1:length(thresh)
          aligns = cellfun(@(x)(find(x>=thresh(i), 1)), all_files(:,2));
          [stack, shift] = stack_images(all_files(:,1), aligns, img_thresh(j));
           avg = nanmean(stack,3);
           residues = bsxfun(@minus, stack, avg);
           residues = residues.^2;
           residues(isnan(residues)) = 0;
           msr(i,j) = sum(residues(:));
        end

        for i=1:3
          [stack, shift] = stack_images(all_files(:,1), timings(:,i), img_thresh(j));
           avg = nanmean(stack,3);
           residues = bsxfun(@minus, stack, avg);
           residues = residues.^2;
           residues(isnan(residues)) = 0;
           msr(i+length(thresh),j) = sum(residues(:));
        end
      end

      figure;imagesc(msr);
      colormap(redgreenmap);
      colorbar
      set(gca, 'YTick', 1:size(msr, 1), ...
               'DataAspectRatio', [1 1 1], ...
               'YTickLabel', {}, ...
               'XTickLabel', {});

      %% + check correlation timings/fraction ??

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
