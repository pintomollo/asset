function find_kymograph(varargin)

  [mymovies, uuid, fitting, opts] = parse_input(varargin{:});

  for i=1:length(mymovies)
    mymovie = mymovies{i};

    if (ischar(mymovie))
      kymos = dir(mymovie);
    elseif (isfield(mymovie, 'mymovie'))
      kymos = mymovie;
    else
      kymos = struct('mymovie', mymovie, 'opts', get_struct('ASSET'));
    end

    for j=1:length(kymos)
      if (isfield(kymos(j), 'mymovie'))
        kymo = kymos(j);
        kymo_name = kymo.mymovie.experiment
      else
        kymo = load(kymos(j).name);
        kymo_name = kymos(j).name;
      end

      try
        if (isfield(kymo, 'mymovie'))
          kymo.opts.recompute = false;
          [fraction, max_width, cell_width, domain, pos] = domain_expansion(kymo.mymovie, kymo.opts);
          times = get_manual_timing(kymo.mymovie, kymo.opts);
          ground_truth = domain(1:times(end), :);

          opts.axes_length = kymo.mymovie.metadata.axes_length_3d;
          opts.reaction_params(end-1,:) = surface2volume(opts.axes_length);
          opts.reaction_params(end, :) = 0.5*ellipse_circum(opts.axes_length);

          kymo_name = kymo.mymovie.experiment;
        elseif (isfield(kymo, 'ground_truth') & isfield(kymo, 'pos'))
          ground_truth = kymo.ground_truth;
          pos = kymo.pos;
          if (isfield(kymo, 'times'))
            times = kymo.times;
          else
            times = NaN(1,3);
          end

          kymo_name = kymos(j).name(1:end-4);
        else
          error('No suitable data for performing the fitting procedure');
        end

        outside = (abs(pos) > opts.reaction_params(end, 1));
        if (any(outside))
          pos = pos(~outside);
          ground_truth = ground_truth(:, ~outside);
        end

        if (fitting.fit_relative)
          fid = fopen(['fitting_rel-' num2str(fitting.fit_flow) '-' num2str(fitting.fit_full) '-' num2str(fitting.parameter_set) '.txt'], 'a');
        else
          fid = fopen(['fitting_adr-' num2str(fitting.fit_flow) '-' num2str(fitting.fit_full) '-' num2str(fitting.parameter_set) '.txt'], 'a');
        end
        fprintf(fid, '%d %s %d\n', fitting.fit_full, kymo_name, fitting.parameter_set);
        fclose(fid);

        boundary = (length(pos)-1)/2;
        fitting.ground_truth = permute([ground_truth(:, 1:boundary+1, :), ground_truth(:,end:-1:boundary+1, :)], [2 1 3]);
        fitting.x_pos = pos(boundary+1:end);
        fitting.t_pos = [0:size(ground_truth, 1)-1]*10;
        fitting.type = kymo_name;

        if (~fitting.fit_full)
          for k=1:size(fitting.ground_truth, 3)
            tmp = fitting.ground_truth(:,:,k);
            tmp = tmp(any(~isnan(tmp), 2), :);
            fitting.ground_truth(end-size(tmp,1)+1:end,:,k) = tmp;
          end

          if (isnan(times(2)))
            if (isnan(times(1)))
              time = 10;
            else
              time = round((times(end) - times(1)) * 0.5);
            end
            fitting.ground_truth = fitting.ground_truth(:, end-time:end, :);
            fitting.t_pos = fitting.t_pos(:, end-time:end);
          else
            fitting.ground_truth = fitting.ground_truth(:, times(2):end, :);
            fitting.t_pos = fitting.t_pos(:, times(2):end);
          end
        end

        if (strncmp(fitting.fitting_type, 'mcmc', 4))
          uuids = metropolis_hastings(fitting, opts);
        elseif (iscell(fitting.fitting_type) | strncmp(fitting.fitting_type, 'all', 3))
          if (iscell(fitting.fitting_type))
            types = fitting.fitting_type;
          else
            types = {'cmaes', 'godlike', 'pso', 'dram'};
          end
          ntypes = length(types);
          types = types(randperm(ntypes));
          uuids = cell(fitting.nfits, 1);
          fitting.nfits = 1;

          for t=1:fitting.nfits
            fitting.fitting_type = types{mod(t-1, ntypes)+1};
            if (strncmp(fitting.fitting_type, 'mcmc', 4))
              uuids{t} = metropolis_hastings(fitting, opts);
            else
              uuids{t} = fit_kymograph(fitting, opts);
            end
          end
        else
          uuids = fit_kymograph(fitting, opts);
        end

        if (fitting.fit_relative)
          fid = fopen(['fitting_rel-' num2str(fitting.fit_flow) '-' num2str(fitting.fit_full) '-' num2str(fitting.parameter_set) '.txt'], 'a');
        else
          fid = fopen(['fitting_adr-' num2str(fitting.fit_flow) '-' num2str(fitting.fit_full) '-' num2str(fitting.parameter_set) '.txt'], 'a');
        end
        for u = 1:length(uuids)
          if (iscell(uuids{u}))
            fprintf(fid, '%s %s OK\n', uuids{u}{1}, kymo_name);
          else
            fprintf(fid, '%s %s OK\n', uuids{u}, kymo_name);
          end
        end
        fclose(fid);

      catch ME

        warning on;
        warning('Error during the fitting procedure');

        if (fitting.fit_relative)
          fid = fopen(['fitting_rel-' num2str(fitting.fit_flow) '-' num2str(fitting.fit_full) '-' num2str(fitting.parameter_set) '.txt'], 'a');
        else
          fid = fopen(['fitting_adr-' num2str(fitting.fit_flow) '-' num2str(fitting.fit_full) '-' num2str(fitting.parameter_set) '.txt'], 'a');
        end
        fprintf(fid, '-1 %s NO\n', kymo_name);
        fclose(fid);

        if (fitting.fit_relative)
          fid = fopen(['fitting_rel-' num2str(fitting.fit_flow) '-' num2str(fitting.fit_full) '-' num2str(fitting.parameter_set) '.log'], 'a');
        else
          fid = fopen(['fitting_adr-' num2str(fitting.fit_flow) '-' num2str(fitting.fit_full) '-' num2str(fitting.parameter_set) '.log'], 'a');
        end
        print_all(fid, ME);
        fclose(fid);

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
      else
        switch varargin{2*i - 1}
          case 'config_modeling'
            conf_mod = varargin{2*i};
          case 'config_fitting'
            conf_fit = varargin{2*i};
          otherwise
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
