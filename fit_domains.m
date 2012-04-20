function fit_domains(fname, incremental, share_work)

  %addpath('./celegans-analysis/');
  %addpath('./celegans-analysis/libraries/');
  %addpath('./pse/');

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
  end

  if (incremental)
    mymovies = struct2cell(mymovies);
    mymovies = mymovies(1,:).';

    fid = fopen('fittings.txt', 'r');
    if (fid > 2)
      dones = {};
      line = fgetl(fid);
      while (ischar(line))
        tmp = regexp(line, ' ', 'split');
        dones(end+1) = {[tmp{2} '.mat']};
        line = fgetl(fid);
      end
      fclose(fid);
      todos = ~ismember(mymovies, dones);
      mymovies = mymovies(todos);
    end

    if (share_work)
      mymovies = mymovies(randi(length(mymovies)));
    end
  end

  for findx = 1:length(mymovies)
    if (incremental)
      kymo = load(mymovies{findx});
    else
      kymo = load(mymovies(findx).name);
    end

    kymo.opts = load_parameters(kymo.opts, 'domain_center.txt');
    domain = imnorm(gather_quantification(kymo.mymovie, kymo.opts));
    kymo.mymovie.data.domain = dynamic_programming(domain, kymo.opts.quantification.params, @weight_symmetry, kymo.opts.quantification.weights, kymo.opts);
    
    nreplicates = 3;
    for t = 0:3
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

      for u = 0:1
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

        for i = 1:nreplicates
          
          options = num2str([findx t u i]); 
          
          if (incremental)
            display([mymovies{findx} ': ' options]);
          else
            display([mymovies(findx).name ': ' options]);
          end
          try
            ml_kymograph(data, kymo);
          catch
            print_all(lasterror)
          end
        end
      end
    end
  end

  return;
end
