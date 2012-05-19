function fit_domains(fname, incremental, share_work)

  addpath('./celegans-analysis/');
  addpath('./celegans-analysis/libraries/');
  addpath('./celegans-analysis/helpers/');
  addpath('./pse/');

  if (nargin < 2)
    incremental = true;
    share_work = true;
  elseif (nargin < 3)
    share_work = true;
  end

  data = pseset('fitting');
  data.fit_noise = 0.075;

  mymovies = [];
  if (ischar(fname))
    mymovies = dir(fname);
  end
  if (isempty(mymovies))
    data.uuid = fname;
    if (ischar(data.uuid))
      data.uuid = str2double(data.uuid);
    end

    mymovies = dir('1056-*_.mat');
  elseif (~islogical(incremental))
    data.uuid = incremental;
    if (ischar(data.uuid))
      data.uuid = str2double(data.uuid);
    end
    incremental = true;
  end

  RandStream.setDefaultStream(RandStream('mt19937ar','Seed',data.uuid));

  replicates = 1:3;
  t_choices = 0:3;
  u_choices = 0:1;

  all_choices = enumerate(t_choices, u_choices, replicates);
  start = 1;

  if (incremental)
    mymovies = struct2cell(mymovies);
    mymovies = mymovies(1,:).';

    fid = fopen('fittings.txt', 'r');
    if (fid > 2)
      ids = {};
      names = {};
      indxs = [];
      running = [];

      line = fgetl(fid);
      while (ischar(line))
        tmp = regexp(line, '^([\d\.]+) (\S+)( \d\s+\d\s+\d\s+\d)?$', 'tokens');
        if (isempty(tmp{1}{3}))
          id = line(1:end-3);
          check = ismember(ids, id);
          running(check) = false;
          indxs(check) = indxs(check) + 1;
        else
          check = ismember(names, tmp{1}{2});
          if (any(check))
            ids(check) = tmp{1}(1);
            running(check) = true;
          else
            ids = [ids; tmp{1}(1)];
            names = [names; tmp{1}(2)];
            indxs = [indxs; 1];
            running = [running; true];;
          end
        end
        line = fgetl(fid);
      end
      fclose(fid);

      names = strcat(names, '.mat');
      finished = (running | indxs >= size(all_choices, 1));
      dones = names(finished);
      names = names(~finished);
      indxs = indxs(~finished);

      todos = ~ismember(mymovies, dones);
      mymovies = mymovies(todos);
    end

    if (share_work)
      mymovies = mymovies(randi(length(mymovies)));

      check = ismember(names, mymovies);
      if (any(check))
        start = indxs(check);
      end
    end
  end

  for findx = 1:length(mymovies)
    if (incremental)
      kymo = load(mymovies{findx});
    else
      kymo = load(mymovies(findx).name);
    end

    if (~isfield(kymo, 'ground_truth'))
      if (all(isnan(get_manual_timing(kymo.mymovie, kymo.opts))))
         error([kymo.mymovie.experiment ' has no valid timing data !']);
      end

      kymo.opts = load_parameters(kymo.opts, 'domain_center.txt');
      domain = imnorm(gather_quantification(kymo.mymovie, kymo.opts));
      kymo.mymovie.data.domain = dynamic_programming(domain, kymo.opts.quantification.params, @weight_symmetry, kymo.opts.quantification.weights, kymo.opts);
    end

    for j = start:size(all_choices, 1)

      t = all_choices(j, 1);
      switch t
        case 0
          data.fit_parameters = [7 12];
          data.lower_bound = zeros(size(data.fit_parameters));
          data.upper_bound = Inf(size(data.fit_parameters));
          data.boundary_sharpness = 35*ones(size(data.fit_parameters));
        case 1
          data.fit_parameters = [7 8 12 13];
          data.lower_bound = zeros(size(data.fit_parameters));
          data.upper_bound = Inf(size(data.fit_parameters));
          data.boundary_sharpness = 35*ones(size(data.fit_parameters));
        case 2
          data.fit_parameters = [7 9 12];
          data.lower_bound = zeros(size(data.fit_parameters));
          data.upper_bound = Inf(size(data.fit_parameters));
          data.boundary_sharpness = 35*ones(size(data.fit_parameters));
        case 3
          data.fit_parameters = [7 8 9 12 13];
          data.lower_bound = zeros(size(data.fit_parameters));
          data.upper_bound = Inf(size(data.fit_parameters));
          data.boundary_sharpness = 35*ones(size(data.fit_parameters));
      end

      u = all_choices(j, 2);
      switch u
        case 0
          % No variability allowed !
        case 1
          data.fit_parameters = [data.fit_parameters 5 10];
          data.lower_bound = [data.lower_bound 0.00688 0.0354];
          data.upper_bound = [data.upper_bound 0.01028 0.0594];
          data.boundary_sharpness = [data.boundary_sharpness 20 20];

          indx = (data.fit_parameters == 9);
          if (any(indx))
            data.lower_bound(indx) = 1.23;
            data.upper_bound(indx) = 1.89;
            data.boundary_sharpness(indx) = 20;
          else
            data.fit_parameters = [data.fit_parameters 9];
            data.lower_bound = [data.lower_bound 1.23];
            data.upper_bound = [data.upper_bound 1.89];
            data.boundary_sharpness = [data.boundary_sharpness 20];
          end
      end

      i = all_choices(j, 3);

      options = num2str([findx t u i]); 

      if (incremental)
        display([mymovies{findx} ': ' options]);
      else
        display([mymovies(findx).name ': ' options]);
      end
      try
        ml_kymograph(data, kymo, options);
      catch
        print_all(lasterror)
      end
    end
  end

  return;
end
