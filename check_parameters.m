function [score, results] = check_parameters(varargin)

  score = Inf;
  results = {};

  [mymovies, uuid, fitting, opts] = parse_input(varargin{:});

  for i=1:length(mymovies)
    mymovie = mymovies{i};

    if (strncmp(fitting.type, 'simulation', 10))
      kymos = struct('ground_truth', -1, 'pos', -1);
    elseif (ischar(mymovie))
      if (strncmp(mymovie(end-2:end), 'txt', 3))
        kymos = textread(mymovie, '%s');
      else
        kymos = dir(mymovie);
        fields = fieldnames(kymos);
        kymos = struct2cell(kymos);
        kymos = kymos(ismember(fields, 'name'), :);
        kymos = kymos(:);
      end
    elseif (isfield(mymovie, 'mymovie'))
      kymos = mymovie;
    else
      kymos = struct('mymovie', mymovie, 'opts', get_struct('ASSET'));
    end

    if (exist(fitting.skip_file, 'file') == 2 && iscell(kymos))
      skip_files = textread(fitting.skip_file, '%s');
      kymos = kymos(~ismember(kymos, skip_files));
    end

    for j=1:length(kymos)
      if (iscell(kymos))
        kymo_name = kymos{j};
        kymo = load(kymo_name);
        kymo_name = kymo_name(1:end-4);
      elseif (isstruct(kymos))
        kymo = kymos(j);
        if (isfield(kymo, 'mymovie'))
          kymo_name = kymo.mymovie.experiment;
        elseif (strncmp(fitting.type, 'simulation', 10))
          kymo_name = 'simulation';
        else
          kymo_name = 'Provided_data';
        end
      else
        error('Unknow object proposed for fitting');
      end

      temps = [];

      try
        if (isfield(kymo, 'mymovie'))
          kymo.opts.recompute = false;
          [fraction, max_width, cell_width, domain, pos] = domain_expansion(kymo.mymovie, kymo.opts);
          times = get_manual_timing(kymo.mymovie, kymo.opts);
          ground_truth = domain(1:times(end), :);

          axes_length = kymo.mymovie.metadata.axes_length_3d;

          kymo_name = kymo.mymovie.experiment;

          temp = regexp(kymo_name, '-(\w*)-', 'tokens');
          temp = str2double(temp{1});

          if (~isfinite(temp))
            fitting.temperature = 23;
          else
            switch temp
              case 14
                fitting.temperature = 20;
              case 3
                fitting.temperature = 13;
              otherwise
                fitting.temperature = temp;
            end
          end
        elseif (isfield(kymo, 'ground_truth') & isfield(kymo, 'pos'))
          ground_truth = kymo.ground_truth;
          pos = kymo.pos;

          if (isfield(kymo, 'times'))
            times = kymo.times;
          else
            times = NaN(1,3);
          end
          if (isfield(kymo, 'temperatures'))
            fitting.temperature = kymo.temperatures;
          else
            temp = regexp(kymo_name, '-(\w*)-', 'tokens');
            temp = str2double(temp{1});

            if (~isfinite(temp))
              fitting.temperature = 23;
            else
              switch temp
                case 14
                  fitting.temperature = 20;
                case 3
                  fitting.temperature = 13;
                otherwise
                  fitting.temperature = temp;
              end
            end
          end
          if (isfield(kymo, 'axes_length'))
            axes_length = kymo.axes_length;
          else
            axes_length = opts.axes_length;
          end
        else
          error('No suitable data for performing the fitting procedure');
        end

        if (~iscell(ground_truth))
          ground_truth = {ground_truth};
          pos = {pos};
        end
        ngroups = length(ground_truth);

        if (size(axes_length, 2) ~= ngroups)
          axes_length = repmat(axes_length(:,1), 1, ngroups);
        end

        if (fitting.fit_relative)
          log_file = ['fitting_rel-' num2str(fitting.fit_flow) '-' num2str(fitting.fit_full) '-' num2str(fitting.parameter_set) '.txt'];
        else
          log_file = ['fitting_adr-' num2str(fitting.fit_flow) '-' num2str(fitting.fit_full) '-' num2str(fitting.parameter_set) '.txt'];
        end

        if (fitting.start_with_best)

          fields = fieldnames(fitting);
          values = struct2cell(fitting);
          filters = {'estimate_n'; 'fit_full'; 'fit_relative'; 'parameter_set';'fit_flow';'fit_model';'aligning_type';'normalize_smooth'};

          good_fields = ismember(fields, filters);

          prev_values = group_ml_results('adr-kymo-*_evol.dat', [{'type', kymo_name} ;[fields(good_fields) values(good_fields)]]);
          if (exist('BestFits', 'dir')==7)
            other_values = group_ml_results(['BestFits' filesep 'adr-kymo-*_evol.dat'], [{'type', kymo_name} ;[fields(good_fields) values(good_fields)]]);
            if (~isempty(other_values))
              if (~isempty(prev_values))
                prev_values{1,2} = [prev_values{1,2}; other_values{1,2}];
              else
                prev_values = other_values;
              end
            end
          end

          if (~isempty(prev_values))
            prev_values = extract_model_parameters(prev_values, true);

            best_score = Inf;
            best_pos = fitting.init_pos;
            best_args = [];

            for k=1:size(prev_values{1,2}, 1)
              if (prev_values{1,2}{k,2}.score < best_score)
                best_score = prev_values{1,2}{k,2}.score;
                best_pos = [prev_values{1,2}{k,2}.params.rate ...
                            prev_values{1,2}{k,2}.params.offset ...
                            prev_values{1,2}{k,2}.params.energy ...
                            prev_values{1,2}{k,2}.params.viscosity ...
                            prev_values{1,2}{k,2}.params.flow ...
                            prev_values{1,2}{k,2}.params.sigma];
              end
            end

            if (~isempty(best_pos))
              fitting.init_pos = best_pos;
            end
          end
        end

        if (~fitting.scale_each_egg)
          [aim_circum, indx] = max(cellfun(@max, pos));
          aim_circum = 2*aim_circum;
          axes_length = repmat(axes_length(:,indx), 1, ngroups);
        end

        axes_length = sort(axes_length, 'descend');
        if (fitting.extrapol_z)
          convert = get_struct('z-correlation');
          axes_length(3,:) = convert.bkg + convert.long_axis*axes_length(1,:) + convert.short_axis*axes_length(2,:);
        end

        for g=1:ngroups

          if (fitting.scale_each_egg)
            aim_circum = 2*max(pos{g});
          end

          if (aim_circum > 0)
            rescale_factor = ellipse_circum(axes_length(:,g), aim_circum, fitting.rescale_length_only);

            if (fitting.rescale_length_only)
              axes_length(1,g) = rescale_factor;
            else
              axes_length(:,g) = axes_length(:,g)*rescale_factor;
            end
          end
          egg_properties(1,g) = surface2volume(axes_length(:,g));
          egg_properties(2,g) = 0.5*ellipse_circum(axes_length(:,g));

          outside = (abs(pos{g}) > egg_properties(2, g));
          if (any(outside))
            pos{g} = pos{g}(~outside);
            ground_truth{g} = ground_truth{g}(:, ~outside, :);
          end

          boundary = (length(pos{g})-1)/2;
          fitting.ground_truth{g} = permute([ground_truth{g}(:, 1:boundary+1, :), ground_truth{g}(:,end:-1:boundary+1, :)], [2 1 3]);
          fitting.x_pos{g} = pos{g}(boundary+1:end);
          fitting.t_pos{g} = [0:size(ground_truth{g}, 1)-1]*10;

          if (~fitting.fit_full)
            for k=1:size(fitting.ground_truth{g}, 3)
              tmp = fitting.ground_truth{g}(:,:,k);
              tmp = tmp(any(~isnan(tmp), 2), :);
              fitting.ground_truth{g}(end-size(tmp,1)+1:end,:,k) = tmp;
            end

            if (isnan(times(2)))
              if (isnan(times(1)))
                time = 10;
              else
                time = round((times(end) - times(1)) * 0.5);
              end
              fitting.ground_truth{g} = fitting.ground_truth{g}(:, end-time:end, :);
              fitting.t_pos{g} = fitting.t_pos{g}(:, end-time:end);
            else
              fitting.ground_truth{g} = fitting.ground_truth{g}(:, times(2):end, :);
              fitting.t_pos{g} = fitting.t_pos{g}(:, times(2):end);
            end
          end
        end
        fitting.egg_properties = egg_properties;
        [score, results] = test_kymograph(fitting, opts);

      catch ME
        print_all(ME)

        continue;
      end
    end
  end

  return;
end

function [mymovies, uuid, fitting, opts] = parse_input(varargin)

  mymovies = '';
  uuid = 1;
  fitting = get_struct('fitting');
  opts = get_struct('modeling');
  opts = load_parameters(opts, 'goehring.txt');

  conf_fit = '';
  conf_mod = '';

  % Check what we got as inputs
  if (nargin > 0)
    mymovies = varargin{1};
    varargin(1) = [];

    if (length(varargin) > 0 & isstruct(varargin{1}))
      if (isfield(varargin{1}, 'reaction_params'))
        opts = varargin{1};
        varargin(1) = [];
      elseif (isfield(varargin{1}, 'parameter_set'))
        fitting = varargin{1};
        varargin(1) = [];
      end
      if (length(varargin) > 0 & isstruct(varargin{1}))
        if (isfield(varargin{1}, 'reaction_params'))
          opts = varargin{1};
          varargin(1) = [];
        elseif (isfield(varargin{1}, 'parameter_set'))
          fitting = varargin{1};
          varargin(1) = [];
        end
      end
    end

    % Now we check that the parameters were provided in pairs
    npairs = length(varargin) / 2;
    if (npairs ~= floor(npairs))
      uuid = str2double(varargin(1));
      varargin(1) = [];
      npairs = floor(npairs);
    end

    % Loop over the pairs of parameters
    for i = 1:npairs
      switch varargin{2*i - 1}
        case 'config_modeling'
          conf_mod = varargin{2*i};
        case 'config_fitting'
          conf_fit = varargin{2*i};
      end
    end

    if (~isempty(conf_mod))
      opts = load_parameters(opts, conf_mod);
    end

    if (~isempty(conf_fit))
      fitting = load_parameters(fitting, conf_fit);
    end

    % Loop over the pairs of parameters
    for i = 1:npairs

      % If the parameter exists in opts we simply assign it the
      % provided value
      if (isfield(fitting, varargin{2*i - 1}))
        fitting.(varargin{2*i - 1}) = varargin{2*i};
      elseif (isfield(opts, varargin{2*i - 1}))
        opts.(varargin{2*i - 1}) = varargin{2*i};
      else
        switch varargin{2*i - 1}
          case 'config_modeling'
            conf_mod = varargin{2*i};
          case 'config_fitting'
            conf_fit = varargin{2*i};
          otherwise
            warning on
            warning(['Property ''' varargin{2*i -1} ''' does not exist. Ignoring']);
        end
      end
    end
  end

  if (~fitting.fit_full)
    opts = load_parameters(opts, 'maintenance.txt');
  end

  if (~iscell(mymovies))
    mymovies = {mymovies};
  end

  return;
end

function [score, results] = test_kymograph(fitting, opts)

  if (nargin == 0)
    fitting = get_struct('fitting');
    opts = get_struct('modeling');
    opts = load_parameters(opts, 'goehring.txt');
    if (~fitting.fit_full)
      opts = load_parameters(opts, 'maintenance.txt');
    end
  elseif (nargin == 1)
    if (isfield(fitting, 'nparticles'))
      opts = fitting;
      fitting = get_struct('fitting');
    else
      opts = get_struct('modeling');
      opts = load_parameters(opts, 'goehring.txt');
      if (~fitting.fit_full)
        opts = load_parameters(opts, 'maintenance.txt');
      end
    end
  else
    orig_opts = get_struct('modeling');
    orig_opts = load_parameters(orig_opts, 'goehring.txt');
    orig_opts = load_parameters(orig_opts, 'custom_model.txt');
    %orig_opts = opts;
    if (fitting.fit_relative)
      orig_opts.reaction_params(3,:) = orig_opts.reaction_params(3,:) .* orig_opts.reaction_params(4,[2 1]);
    end
    orig_params = [orig_opts.diffusion_params; ...
                orig_opts.reaction_params];

    orig_scaling = 10.^(floor(log10(orig_params)));
    orig_scaling(orig_params == 0) = 1;
  end

  if (fitting.fit_relative)
    opts.reaction_params(3,:) = opts.reaction_params(3,:) ./ opts.reaction_params(4,[2 1]);
  end

  if (~isempty(fitting.simulation_parameters))
    fitting.simulation_parameters = reshape(fitting.simulation_parameters, [], 2);

    opts.diffusion_params = fitting.simulation_parameters(1,:);
    opts.reaction_params = fitting.simulation_parameters(2:end,:);
  end

  if (~isempty(fitting.flow_size) && size(opts.advection_params, 2) ~= fitting.flow_size(end))
    tmp_opts = get_struct('modeling');
    tmp_opts = load_parameters(tmp_opts, 'goehring.txt');

    opts.advection_params = tmp_opts.advection_params;
    opts.user_data = tmp_opts.user_data;
    opts.flow_temperature = tmp_opts.flow_temperature;
    opts.tmax = tmp_opts.tmax;
  end

  ngroups = length(fitting.ground_truth);
  results = cell(ngroups, 2);

  all_params = false;
  almost_all_params = false;

  noffsets = ngroups*strncmp(fitting.aligning_type, 'fitting', 7);
  [temperatures, junk, temp_indx] = unique(fitting.temperature);
  good_visc = (temperatures ~= opts.reaction_temperature);
  bad_flow = (temperatures == opts.flow_temperature);
  nvisc = sum(good_visc);
  ntemps = length(temperatures);
  viscosities = ones(1, ntemps);
  has_fixed_parameters = false;

  ml_params = [opts.diffusion_params; ...
                opts.reaction_params];

  ml_params(end-1,:) = fitting.egg_properties(1,1);
  ml_params(end, :) = fitting.egg_properties(2,1);
  half = ml_params(end, 1);

  flow_scaling = 1;
  if (fitting.scale_flow || ~isempty(opts.scale_params))
    data_flow_scaling = false;
    if (numel(opts.scale_params)>1 && opts.scale_params(2))
      tmp_opts = get_struct('modeling');
      tmp_opts = load_parameters(opts, 'test_data_flow');
      tmp_conv = get_struct('z-correlation');
      tmp_factor = surface2volume([tmp_opts.axes_length(1:2); tmp_conv.bkg + tmp_conv.long_axis*tmp_opts.axes_length(1) + tmp_conv.short_axis*tmp_opts.axes_length(2)]);
      flow_scale_factor = (tmp_factor ./ fitting.egg_properties(1,:)) - 1;
      data_flow_scaling = true;
    else
      tmp_opts = get_struct('modeling');
      tmp_conv = get_struct('conversion');
      tmp_factor = surface2volume(tmp_opts.axes_length .* [tmp_conv.maintenance;1;1]);
      %flow_scale_factor = (fitting.egg_properties(1,:) / tmp_factor) - 1;
      flow_scale_factor = (tmp_factor ./ fitting.egg_properties(1,:)) - 1;
    end
  end

  [fit_params, fit_energy, fit_temperatures, fit_viscosity] = model_description(fitting.parameter_set);

  % In these case, we always need a model
  if (mod(fitting.parameter_set, 10) < 3)
    fitting.fit_model = true;
  end

  nrates = length(fit_params);
  nenergy = length(fit_energy);

  if (~fitting.fit_model && fit_temperatures)
    fit_energy = ones(1, nenergy*nvisc);
    fitting.fit_flow = any(~good_visc);
  end

  if (~isempty(fitting.init_pos))
    fitting.init_pos = fitting.init_pos(:).';

    if (iscell(fitting.init_pos))
      tmp_fit = fitting;

      for c=1:length(fitting.init_pos)
        tmp_fit.init_pos = fitting.init_pos{c};
        score(c) = test_kymograph(tmp_fit, opts);
      end

      return;
    elseif ((nrates + numel(fit_energy) + noffsets + fitting.fit_flow + fit_viscosity*nvisc + fitting.scale_flow) == numel(fitting.init_pos))

      if (nrates > 0)
        fitting.init_pos(1:nrates) = fitting.init_pos(1:nrates) ./ orig_scaling(fit_params);
      end
      fitting.init_pos(nrates+1:nrates+noffsets) = fitting.init_pos(nrates+1:nrates+noffsets) ./ fitting.offset_scaling;

      all_params = true;
    elseif ((nrates + numel(fit_energy) + fitting.fit_flow + fit_viscosity*nvisc + fitting.scale_flow) == numel(fitting.init_pos))

      if (nrates > 0)
        fitting.init_pos(1:nrates) = fitting.init_pos(1:nrates) ./ orig_scaling(fit_params);
      end
      almost_all_params = true;
    elseif ((nrates + numel(fit_energy) + noffsets + fit_viscosity*nvisc) == numel(fitting.init_pos))

      if (nrates > 0)
        fitting.init_pos(1:nrates) = fitting.init_pos(1:nrates) ./ orig_scaling(fit_params);
      end

      fitting.init_pos(nrates+1:nrates+noffsets) = fitting.init_pos(nrates+1:nrates+noffsets) ./ fitting.offset_scaling;
    elseif (nrates == numel(fitting.init_pos))
      fitting.init_pos = fitting.init_pos ./ orig_scaling(fit_params);
      fitting.init_pos = [fitting.init_pos fit_energy];
    elseif (4 + noffsets == numel(fitting.init_pos))
      warning('Assuming only the 4 main ones were provided !')
      tmp_params = ml_params;
      tmp_params([4 5 12 13]) = fitting.init_pos(1:4);
      tmp_params = tmp_params ./ orig_scaling;

      fitting.init_pos = [tmp_params(fit_params) fitting.init_pos(5:end)/fitting.offset_scaling fit_energy];
    elseif (4 == numel(fitting.init_pos))
      warning('Assuming only the 4 main ones were provided !')
      tmp_params = ml_params;
      tmp_params([4 5 12 13]) = fitting.init_pos(1:4);
      tmp_params = tmp_params ./ orig_scaling;

      fitting.init_pos = [tmp_params(fit_params) fit_energy];
    else
      warning('The provided initial position does not correspond to the dimensionality of the fit, ignoring it.');
      keyboard
      fitting.init_pos = [];
    end
  end

  if (~isempty(fitting.init_pos) && fitting.fit_relative)
    indx1 = find(fit_params == 4, 1);
    if (~isempty(indx1))
      indx2 = find(fit_params == 13, 1);
      if (~isempty(indx2))
        fitting.init_pos(indx1) = fitting.init_pos(indx1) / fitting.init_pos(indx2);
      else
        fitting.init_pos(indx1) = fitting.init_pos(indx1) / opts.reaction_params(4,2);
      end
    end

    indx1 = find(fit_params == 12, 1);
    if (~isempty(indx1))
      indx2 = find(fit_params == 5, 1);
      if (~isempty(indx2))
        fitting.init_pos(indx1) = fitting.init_pos(indx1) / fitting.init_pos(indx2);
      else
        fitting.init_pos(indx1) = fitting.init_pos(indx1) / opts.reaction_params(4,1);
      end
    end
  end

  if (exist('rng'))
    rng(now + cputime, 'twister');
  else
    RandStream.setDefaultStream(RandStream('mt19937ar','seed',now + cputime));
  end
  normalization_done = false;

  if (strncmp(fitting.fitting_type, 'dram', 4))
    drscale  = 2; 
    adaptint = 0;
  end

  if (~strncmp(fitting.type, 'simulation', 10))
    if (isempty(fitting.ground_truth))
      error('Data is missing for the fitting !');
    end

    opts.boundaries = [0 fitting.egg_properties(2,1)];
    opts.x_step = diff(opts.boundaries)/(opts.nparticles-1);
  end

  flow = opts.advection_params;
  if (size(flow, 1) ~= opts.nparticles)
    [X, Y] = meshgrid([1:size(flow, 2)], 1+([0:opts.nparticles-1]*(size(flow, 1)-1)/(opts.nparticles-1)).');
    flow = bilinear_mex(flow, X, Y, [2 2]);
  end

  for g=1:ngroups

    opts.reaction_params(end-1,:) = fitting.egg_properties(1,g);
    opts.reaction_params(end, :) = fitting.egg_properties(2,g);

    opts.boundaries = [0 fitting.egg_properties(2,g)];
    opts.x_step = diff(opts.boundaries)/(opts.nparticles-1);

    is_data = ~(strncmp(fitting.type, 'simulation', 10));
    if (~is_data)

        if (~isempty(opts.init_func) & isa(opts.init_func, 'function_handle'))
          x0 = opts.init_func(opts, fitting.fit_relative);
        else
          x0 = opts.init_params;
          x0 = repmat(x0, [opts.nparticles, 1]);
        end

        if (fitting.fit_relative)
           [fitting.ground_truth, fitting.t_pos] = simulate_model_real(x0, ml_params, opts.x_step, opts.tmax*0.75, opts.time_step, opts.output_rate, flow, opts.user_data, opts.max_iter);
        else
           [fitting.ground_truth, fitting.t_pos] = simulate_model_mix(x0, ml_params, opts.x_step, opts.tmax*0.75, opts.time_step, opts.output_rate, flow, opts.user_data, opts.max_iter);
        end

      fitting.ground_truth = fitting.ground_truth((end/2)+1:end, :);
      fitting.ground_truth = {[fitting.ground_truth; fitting.ground_truth]};
      fitting.x_pos = {[0:opts.nparticles-1] * opts.x_step};

    elseif (~isempty(fitting.t_pos))
      if (numel(fitting.t_pos{g}) == 1)
        opts.output_rate(g) = fitting.t_pos{g};
      else
        opts.output_rate(g) = median(diff(fitting.t_pos{g}));
      end
    end

    if (~fitting.fit_full)
      nmaintenance = min(size(fitting.ground_truth{g}, 2), 50) - 1;
      fitting.ground_truth{g} = fitting.ground_truth{g}(:, end-nmaintenance:end, :);
      fitting.t_pos{g} = fitting.t_pos{g}(end-nmaintenance:end);
      fitting.aligning_type = 'end';
      fitting.fit_flow = false;
    end

    [n,m,o] = size(fitting.ground_truth{g});
    size_data(g,1:3) = [n,m,o];

    opts_expansion = load_parameters(get_struct('ASSET'), 'domain_expansion.txt');
    [f, frac_width, full_width] = domain_expansion(mymean(fitting.ground_truth{g}(1:end/2, :, :), 3).', mymean(fitting.ground_truth{g}((end/2)+1:end, :, :), 3).', size_data(g,1)/2, size_data(g,2), opts_expansion);

    [junk, relative_fraction{g}] = synchronize_domains(f, f);
    tmp_gfrac = f*frac_width;
    tmp_gfrac(isnan(tmp_gfrac)) = 0;
    gfraction{g} = tmp_gfrac;

    if (strncmp(fitting.scale_type, 'normalize', 10))
      fitting.ground_truth{g} = normalize_domain(fitting.ground_truth{g}, f*frac_width, opts_expansion, is_data, fitting.normalize_smooth);

      normalization_done = true;
    end

    if (strncmp(fitting.aligning_type, 'domain_width', 12))
      frac_indx(g) = find(f > fitting.fraction, 1, 'first');
    end

    if (~normalization_done & strncmp(fitting.scale_type, 'normalize', 10))
      opts_expansion = load_parameters(get_struct('ASSET'), 'domain_expansion.txt');
      [f, frac_width, full_width] = domain_expansion(mymean(fitting.ground_truth{g}(1:end/2, :, :), 3).', mymean(fitting.ground_truth{g}((end/2)+1:end, :, :), 3).', size_data(g,1)/2, size_data(g,2), opts_expansion);
      if (~fitting.fit_full)
        f(:) = 1;
      end
      fitting.ground_truth{g} = normalize_domain(fitting.ground_truth{g}, f*frac_width, opts_expansion, is_data, fitting.normalize_smooth);
      normalization_done = true;

      tmp_gfrac = f*frac_width;
      tmp_gfrac(isnan(tmp_gfrac)) = 0;
      gfraction{g} = tmp_gfrac;
    end

    if (strncmp(fitting.aligning_type, 'lsr', 3))
      sindx(g) = round(size(fitting.ground_truth{g}, 2) / 2);
      mean_ground_truth{g} = nanmean(fitting.ground_truth{g}, 3);
    end

    multi_data(g) = (size_data(g,3) > 1);
    if (multi_data(g))
      nlayers(g) = size_data(g,3);
    end

    if (~fitting.integrate_sigma)
      noise = estimate_noise(fitting.ground_truth{g});
      estim_sigma(g) = mymean(noise(:,2));
    end

    linear_truth{g} = fitting.ground_truth{g}(:);
    linear_goods{g} = isfinite(linear_truth{g});
    simul_pos = [0:opts.nparticles-1].';
    penalty(g) = ((max(linear_truth{g}(linear_goods{g})))^2)*opts.nparticles;
    ndata(g) = length(fitting.x_pos{g});

    if (fit_temperatures || ~isempty(opts.temperature_params))

      kB = 8.6173324e-5;
      C2K = 273.15;

      switch length(opts.temperature_params)
        case 0
          E = ones(3,2)*0.65;
          flow_E = 0.65;
          visc_E = 0.65;
        case 1
          E = ones(3,2)*opts.temperature_params;
          flow_E = E(1);
          visc_E = E(1);
        case 2
          E = ones(3,2)*opts.temperature_params(1);
          flow_E = opts.temperature_params(1);
          visc_E = opts.temperature_params(2);
        case 3
          E = ones(3,2)*opts.temperature_params(1);
          flow_E = opts.temperature_params(2);
          visc_E = opts.temperature_params(3);
        otherwise
          if (fit_temperatures)
            warning('Temperature parameters provided in the configuration file does not correspond to a temeprature model, ignoring them');

            E = ones(3,2)*0.65;
            flow_E = 0.65;
            visc_E = 0.65;
          else
            opts.temperature_params = [];
            E = NaN;
          end
      end

      if (numel(opts.temperature_params) == numel(fit_energy))
        fit_energy = opts.temperature_params;
      end

      if (isfinite(E))
        fit_temperatures = true;

        diff_ratio = ((fitting.temperature(g)+C2K) / (opts.reaction_temperature+C2K)) * ...
                     exp(-(visc_E/kB)*((1/(fitting.temperature(g)+C2K)) - (1/(opts.reaction_temperature+C2K))));
        ratio = exp(-(E/kB)*((1/(fitting.temperature(g)+C2K)) - (1/(opts.reaction_temperature+C2K))));
        flow_ratio = exp(-(flow_E/kB)*((1/(fitting.temperature(g)+C2K)) - (1/(opts.flow_temperature+C2K))));

        temp_scale{g} = [ones(1,2)*diff_ratio; ratio; ones(4,2)];
        flow_scale(g) = flow_ratio;
      end
    else
      E = [];
    end
  end

  rescaling = orig_scaling;
  ml_params = ml_params ./ rescaling;

  if (fitting.scale_each_egg)
    egg_properties = bsxfun(@rdivide, fitting.egg_properties, rescaling(end-1:end,1));
  end

  if (strncmp(fitting.type, 'simulation', 10))
    noiseless = fitting.ground_truth;
    for g=1:ngroups
      range_data(g) = fitting.data_noise * range(linear_truth{g});
    end
  end

  warning off;
  useful_data = [];

  if (fitting.cross_measure & multi_data)
    fitting.nfits = 2*fitting.nfits;
    orig_fit = fitting;
  end

    if (fitting.cross_measure)
      fitting = orig_fit;
      for g=1:ngroups
        if (mod(f, 2) == 1)
          selection = randperm(size(fitting.ground_truth{g}, 3));
          nlayers(g) = round(length(selection)/2);
          fitting.ground_truth{g} = fitting.ground_truth{g}(:,:,selection(1:nlayers(g)));
        else
          fitting.ground_truth{g} = fitting.ground_truth{g}(:,:,selection(nlayers(g)+1:end));
          nlayers(g) = length(selection) - nlayers(g);
        end
        linear_truth{g} = fitting.ground_truth(:);
        linear_goods = isfinite(linear_truth{g});
      end
    end

    tmp_params = ml_params;

    if (isempty(fitting.init_pos))
      p0 = ml_params(fit_params);
    else
      p0 = fitting.init_pos;
    end

    if (~fit_temperatures)
      flow_scale = 1;
    elseif (isempty(fitting.init_pos))
      p0 = [p0 fit_energy];
    end

    if (~fit_viscosity)
      viscosities(:) = 1;
    elseif ~(all_params || almost_all_params)
      p0 = [p0 viscosities(good_visc)];
    end

    if (fitting.fit_flow & ~(all_params || almost_all_params))
      p0 = [p0 1];
    else
      curr_flow_scale = 1;
    end

    if (fitting.fit_sigma)
      p0 = [p0 estim_sigma];
    end

    if (fitting.scale_flow & ~(all_params || almost_all_params))
      if (~isempty(opts.scale_params))
        p0 = [p0 opts.scale_params(1)];
      else
        p0 = [p0 1];
      end
    end

    correct = 0;
    tmp_params = ml_params;
    tmp_opts = opts;

    tmp_opts.diffusion_params = tmp_params(1, :) .* rescaling(1, :);
    tmp_opts.reaction_params = tmp_params(2:end, :) .* rescaling(2:end, :);

    [x0] = opts.init_func(tmp_opts, fitting.fit_relative);

    ndecimals = -min(floor(real(log10(fitting.tolerance/10))), 0);
    if (~isfinite(ndecimals))
      ndecimals = 0;
    end

    each_full_error = NaN(1, ngroups);
    for g=1:ngroups
      if (fitting.estimate_n)
        ns = identify_n(fitting.ground_truth{g});
        nobs(g) = sum([ns(:,1); ns(:,3)]);
        norm_coeff(g) = sum(linear_goods{g})/nobs(g);
      else
        nobs(g, 1:2) = [sum(linear_goods{g}) size_data(g,2)];
        norm_coeff(g) = 1;
      end

      curr_score_weights(g) = fitting.score_weights*length(gfraction{g})/nobs(g,1);

      if (fitting.integrate_sigma)
        if (fitting.pixels_only)
          each_full_error(g) = (0.5*nobs(g,1))*log(penalty(g) * size_data(g,2) * 10);
        else
          each_full_error(g) = (0.5*nobs(g,2))*log(curr_score_weights(g)*penalty(g)*size_data(g,2)*10 + size_data(g,2)*10);
        end
      else
        each_full_error(g) = (curr_score_weights(g) * penalty(g) * size_data(g,2) * 10 + prod(size_data(g,1:2))) / (2*estim_sigma(g).^2);
      end

      if (strncmp(fitting.type, 'simulation', 10))
        fitting.ground_truth{g} = noiseless{g} + range_data(g)*randn(size_data(g,:));
      end
    end
    full_error = sum(each_full_error);

    if (strncmp(fitting.aligning_type, 'fitting', 7) && length(p0) <= (nrates + numel(fit_energy) + fitting.fit_flow + nvisc*fit_viscosity + fitting.scale_flow))
      nparams = length(p0);
      fitting.aligning_type = 'domain';
      [junk, offsets] = error_function(p0(:));
      fitting.aligning_type = 'fitting';
      p0 = [p0(1:nrates) offsets/fitting.offset_scaling p0(nrates+1:end)];
    end

    if (~fitting.scale_flow && ~isempty(opts.scale_params))
      p0 = [p0 opts.scale_params(1)];
      fitting.scale_flow = true;
      fitting.fixed_parameter = [fitting.fixed_parameter length(p0)];
    end

    nparams = length(p0);

    tmp_fit = fitting;
    tmp_fit.ground_truth = [];
    tmp_fit.flow_size = size(flow);
    tmp_fit.rescale_factor = rescaling(fit_params);
    tmp_fit.simulation_parameters = ml_params .* rescaling;

    fitting.scale_each_egg = (fitting.scale_each_egg && (ngroups > 1));

    score = error_function(p0(:));

    warning on;

  return;

  function [err_all, offsets] = error_function(varargin)

    p_all = varargin{1};
    p_all(~fitting.fixed_parameter) = roundn(p_all(~fitting.fixed_parameter), -ndecimals);

    [curr_nparams, nevals] = size(p_all);
    flip = false;
    if (nevals == nparams)
      flip = true;
      p_all = p_all.';
      nevals = curr_nparams; 
    end
    err_all = NaN(ngroups, nevals);

    offsets = zeros(1, ngroups);

    for i = 1:nevals
      correct = 3;

      for g=1:ngroups
        curr_p = p_all(:,i);

        tmp_params = ml_params;

        if (fitting.scale_each_egg)
          tmp_params(end-1,:) = egg_properties(1, g);
          tmp_params(end, :) = egg_properties(2, g);

          opts.boundaries = [0 fitting.egg_properties(2,g)];
          opts.x_step = diff(opts.boundaries)/(opts.nparticles-1);

          half = opts.boundaries(2);
        end

        tmp_params(fit_params) = abs(curr_p(1:nrates));
        more_params = curr_p(nrates+1:end);

        if (fitting.scale_flow)
          if (data_flow_scaling)
            flow_scaling = (0.34*flow_scale_factor(g) + 1);
          else
            flow_scaling = (-abs(more_params(end))*flow_scale_factor(g) + 1);
          end
          more_params = more_params(1:end-1);
        end

        if (fitting.fit_sigma)
          estim_sigma = abs(more_params(end));
          more_params = more_params(1:end-1);
        end

        if (fitting.fit_flow)
          curr_flow_scale = abs(more_params(end));
          more_params = more_params(1:end-1);
        else
          curr_flow_scale = 1;
        end
        curr_flow_scale = curr_flow_scale*flow_scaling;

        if (fit_viscosity)
          curr_visc = ones(1,ntemps);
          curr_visc(good_visc) = abs(more_params(end-nvisc+1:end));
          more_params = more_params(1:end-nvisc);
        else
          curr_visc = viscosities;
        end

        if (fit_temperatures)
          if (isempty(fit_energy))
            tmp_params = tmp_params .* temp_scale{g};
            curr_flow_scale = curr_flow_scale * flow_scale(g);
          else
            if (fitting.fit_model)
              switch nenergy
                case 1
                  E = ones(3,2)*abs(more_params(end));
                  flow_E = E(1);
                  more_params = more_params(1:end-1);
                case 2
                  E = ones(3,2)*abs(more_params(end-1));
                  flow_E = abs(more_params(end));
                  more_params = more_params(1:end-2);
                case 3
                  E = ones(3,2)*abs(more_params(end-2));
                  flow_E = abs(more_params(end-1));
                  visc_E = abs(more_params(end));
                  more_params = more_params(1:end-3);
                case 4
                  E = abs(more_params(end-3:end-1));
                  E = E(:);
                  E = E(:,[1 1]);
                  flow_E = abs(more_params(end));
                  more_params = more_params(1:end-4);
                case 7
                  E = abs(more_params(end-6:end-1));
                  E = E(:);
                  E = [E(1:3,1) E(4:6,1)];
                  flow_E = abs(more_params(end));
                  more_params = more_params(1:end-7);
                otherwise
                  error(['No temperature model with ' num2str(nenergy) ' parameters has been implemented.']);
              end

              diff_ratio = ((fitting.temperature(g)+C2K) / (opts.reaction_temperature+C2K)) * ...
                           exp(-(visc_E/kB)*((1/(fitting.temperature(g)+C2K)) - (1/(opts.reaction_temperature+C2K))));
              ratio = exp(-(E/kB)*((1/(fitting.temperature(g)+C2K)) - (1/(opts.reaction_temperature+C2K))));
              flow_ratio = exp(-(flow_E/kB)*((1/(fitting.temperature(g)+C2K)) - (1/(opts.flow_temperature+C2K))));
            else

              diff_ratio = 1;

              curr_ratios = ones(nenergy,ntemps);
              curr_ratios(:,good_visc) = reshape(abs(more_params(end-nenergy*nvisc+1:end)), nenergy, nvisc);
              more_params = more_params(1:end-nenergy*nvisc);

              eindx = temp_indx(g);

              switch nenergy
                case 1
                  ratio = ones(3,2)*curr_ratios(1,eindx);
                  flow_ratio = ratio(1);
                case 2
                  ratio = ones(3,2)*curr_ratios(1,eindx);
                  flow_ratio = curr_ratios(2, eindx);
                case 4
                  ratio = curr_ratios(1:3, [eindx eindx]);
                  flow_ratio = curr_ratios(4, eindx);
                case 7
                  ratio = [curr_ratios(1:3, eindx) curr_ratios(4:6, eindx)];
                  flow_ratio = curr_ratios(7, eindx);
                otherwise
                  error(['No temperature scaling with ' num2str(nenergy) ' parameters has been implemented.']);
              end

              if (bad_flow(eindx))
                flow_ratio = 1;
                curr_flow_scale = 1;
              elseif (good_visc(eindx))
                curr_flow_scale = 1;
              end
            end

            tmp_params = tmp_params .* [ones(1,2)*diff_ratio; ratio; ones(4,2)];
            curr_flow_scale = curr_flow_scale * flow_ratio;
          end
        end
        if (fit_viscosity)
          tmp_params(1,:) = tmp_params(1,:)*curr_visc(temp_indx(g));
        end

        if (opts.restart_init)
          opts.diffusion_params = tmp_params(1, :) .* rescaling(1, :);
          opts.reaction_params = tmp_params(2:end, :) .* rescaling(2:end, :);
          [x0, correct] = opts.init_func(opts, fitting.fit_relative);

          if (correct == 0)
            err_all(g,i) = Inf;
            continue;
          end
        end

        normalization_done = false;
        if (fitting.fit_relative)
          [res, t] = simulate_model_real(x0, tmp_params .* rescaling, opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow * curr_flow_scale, opts.user_data, opts.max_iter);
        else
          [res, t] = simulate_model_mix(x0, tmp_params .* rescaling, opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow * curr_flow_scale, opts.user_data, opts.max_iter);
        end

        res = res((end/2)+1:end, :);
        if (opts.nparticles ~= ndata)
          res = flipud(interp1q(simul_pos*opts.x_step, flipud(res), fitting.x_pos{g}.'));
        end

        [f, fwidth] = domain_expansion(res.', size(res, 1), size(res,2), opts_expansion);
        fraction = f*fwidth;
        res = [res; res];

        switch fitting.aligning_type
          case 'best'
            tmp_frac = fraction;
            tmp_frac(isnan(fraction)) = 0;

            tmp_res = res;
            tmp_res(isnan(res)) = 0;

            tmp_truth = inpaint_nans(mymean(fitting.ground_truth{g}, 3));

            if (length(t) < size_data(g,2))
              cc = normxcorr2(tmp_res, tmp_truth); 
              cc2 = normxcorr2(tmp_frac, gfraction{g}); 

              [max_cc, imax] = max(cc(size(res, 1), 1:2*size(res,2)));
              [max_cc, imax2] = max(cc2(1:2*size(res,2)));
              corr_offset = ((imax-size(res,2))*fitting.score_weights + (imax2-size(res,2))) / (1 + fitting.score_weights);
            else
              cc = normxcorr2(tmp_truth, tmp_res); 
              cc2 = normxcorr2(gfraction{g}, tmp_frac); 

              [max_cc, imax] = max(cc(size_data(g,1), 1:2*size_data(g,2)));
              [max_cc, imax2] = max(cc2(1:2*size_data(g,2)));
              corr_offset = ((imax-size_data(g,2))*fitting.score_weights + (imax2-size_data(g,2))) / (1 + fitting.score_weights);
            end
          case 'domain'
            if (size(res,2) <= 10)
              err_all(g,i) = Inf;
              continue;
            else

              if (strncmp(fitting.scale_type, 'normalize', 10))
                res = normalize_domain(res, fraction, opts_expansion, false, fitting.normalize_smooth);
                normalization_done = true;
              end
            end

            if (isnan(f(end)))
              err_all(g,i) = Inf;
              continue;
            end

            corr_offset = synchronize_domains(relative_fraction{g}, f);
            disp_args = corr_offset;
            corr_offset = corr_offset(1);

            if (isempty(corr_offset))
              corr_offset = NaN;
            end

            offsets(g) = corr_offset;
          case 'domain_width'
            if (size(res,2) <= 10)
              err_all(g,i) = Inf;
              continue;
            else

              if (strncmp(fitting.scale_type, 'normalize', 10))
                res = normalize_domain(res, fraction, opts_expansion, false, fitting.normalize_smooth);
                normalization_done = true;
              end
            end

            if (isnan(f(end)))
              err_all(g,i) = Inf;
              continue;
            end
            findx = find(f > fitting.fraction, 1, 'first');
            corr_offset = findx - frac_indx(g);
          case 'fitting'
            corr_offset = more_params(g)*fitting.offset_scaling;
          case 'end'
            corr_offset = size_data(g,2) - size(res, 2) + 1;
          case 'lsr'
            rindx = min(sindx(g), size(res, 2));

            [junk, junk2, rindx] = find_min_residue(mean_ground_truth{g}, sindx(g), res, rindx, 0.95);
            corr_offset = sindx(g) - rindx;
        end

        if (size_data(g,2) - corr_offset < 10)
          err_all(g,i) = each_full_error(g);
        else
          if (~normalization_done & strncmp(fitting.scale_type, 'normalize', 10))
            if (~fitting.fit_full)
              fraction(:) = fwidth;
            end
            res = normalize_domain(res, fraction, opts_expansion, false, fitting.normalize_smooth);
            normalization_done = true;
          end

          fraction(isnan(fraction)) = 0;
          res(isnan(res)) = 0;

          res = interp2([res; fraction.'], [1:length(gfraction{g})].'+corr_offset, [1:size(fitting.ground_truth{g},1)+1]);
          fraction = res(end,:).';
          res = res(1:end-1,:);

          fraction(isnan(fraction)) = 0;
          res(isnan(res)) = 0;

          results{g,1} = res;
          results{g,2} = fraction;

          if (multi_data)
            res = repmat(res, [1 1 nlayers(g)]);
          end

          if (~normalization_done & strncmp(fitting.scale_type, 'best', 4))
            gres = ~isnan(res(:)) & linear_goods{g};
            c = [ones(sum(gres), 1), res(gres)] \ linear_truth{g}(gres);

            if (c(2) <= 0)
              err_all(g,i) = Inf;
              continue;
            end

            res = c(1) + c(2)*res;
          end

          tmp_err = (fitting.ground_truth{g} - res).^2 / norm_coeff(g);

          tmp_frac = 1./(1+exp(-15*(((gfraction{g} - fraction)/half).^2-0.5)));

          if (fitting.integrate_sigma)
            % Integrated the sigma out from the gaussian error function
            if (fitting.pixels_only)
              err_all(g,i) = (0.5*nobs(g,1))*log(sum(tmp_err(linear_goods{g})));
            else
              err_all(g,i) = (0.5*nobs(g,2))*log(curr_score_weights(g)*sum(tmp_err(linear_goods{g})) + sum(tmp_frac));
            end
          else
            if (fitting.fit_sigma)
              err_all(g,i) = (curr_score_weights(g)*sum(tmp_err(linear_goods{g})) + sum(tmp_frac)) / (2*estim_sigma(g).^2) + nobs(g)*log(estim_sigma) ;
            else
              err_all(g,i) = (curr_score_weights(g)*sum(tmp_err(linear_goods{g})) + sum(tmp_frac)) / (2*estim_sigma(g).^2);
            end
          end
        end
      end
    end

    err_all(isnan(err_all)) = Inf;
    err_all = sum(err_all, 1);

    if (flip)
      err_all = err_all.';
    end

    if (all(isinf(err_all)))
      % ML crashes when all values are non-numerical
      err_all(1) = full_error*sign(err_all(1));
    end

    return;
  end
end

function domain = normalize_domain(domain, path, opts, has_noise, do_min_max)

  prct_thresh = 5;
  path = path/opts.quantification.resolution;
  [h, w, f] = size(domain);
  h = h/2;
  pos_mat = repmat([1:h].', 1, w);
  mask = bsxfun(@le, pos_mat, path.');
  mask = repmat(flipud(mask), [2, 1]);

  nplanes = size(domain, 3);

  if (do_min_max && has_noise)
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

  return;
end
