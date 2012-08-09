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
      else
        kymo = load(kymos(j).name);
      end

      try
        fid = fopen(['fitting_adr-' num2str(fitting.fit_full) '-' sum2str(fitting.parameter_set) '.txt'], 'a');
        fprintf(fid, '%d %s %d\n', fitting.fit_full, kymo.mymovie.experiment, fitting.parameter_set);
        fclose(fid);

        domain = imnorm(gather_quantification(kymo.mymovie, kymo.opts));
        kymo.opts = load_parameters(kymo.opts, 'domain_center.txt');
        kymo.mymovie.data.domain = dynamic_programming(domain, kymo.opts.quantification.params, @weight_symmetry, kymo.opts.quantification.weights, kymo.opts);

        [ground_truth, junk, pos, indx] = align_domain(kymo.mymovie, kymo.opts);
        valids = any(isnan(ground_truth), 1);

        first = indx - find(valids(1:indx), 1, 'last') - 1;
        last = find(valids(indx:end), 1, 'first') - 2;
        boundary = min(first, last);

        times = get_manual_timing(kymo.mymovie, kymo.opts);

        ground_truth = ground_truth(:,[-boundary:boundary]+indx);
        ground_truth = ground_truth(1:times(end), :);
        g = [ground_truth(:,1:boundary+1), fliplr(ground_truth(:,boundary+1:end))].';
        p = pos(indx:indx+boundary);

        fitting.ground_truth = g;
        fitting.x_pos = p;
        fitting.t_pos = [0:size(g, 2)-1]*10;
        fitting.type = 'data';

        if (~fitting.fit_full)
          fitting.ground_truth = fitting.ground_truth(:, times(2):end);
          fitting.t_pos = fitting.t_pos(:, times(2):end);
        end

        uuids = fit_kymograph(fitting, opts);

        fid = fopen(['fitting_adr-' num2str(fitting.fit_full) '-' sum2str(fitting.parameter_set) '.txt'], 'a');
        for u = 1:length(uuids)
          fprintf(fid, '%s %s OK\n', uuids{u}, kymo.mymovie.experiment);
        end
        fclose(fid);

      catch ME
        print_all(ME);

        fid = fopen(['fitting_adr-' num2str(fitting.fit_full) '-' sum2str(fitting.parameter_set) '.txt'], 'a');
        fprintf(fid, '-1 %s NO\n', kymo.mymovie.experiment);
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

    % Now we check that the parameters were provided in pairs
    npairs = length(varargin) / 2;
    if (npairs ~= floor(npairs))
      uuid = str2double(varargin(1));
      varargin(1) = [];
      npairs = floor(npairs);
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

  if (~isempty(conf_mod))
    opts = load_parameters(opts, conf_mod);
  end

  if (~isempty(conf_fit))
    fitting = load_parameters(fitting, conf_fit);
  end

  if (~iscell(mymovies))
    mymovies = {mymovies};
  end

  return;
end
