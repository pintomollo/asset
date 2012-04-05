function fit_domains(fname)

  mymovies = dir(fname);
  data = pseset('fitting');
  data.fit_noise = 0.075;

  for indx = 1:length(mymovies)
    kymo = load(mymovies(indx).name);
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
          display([mymovies(indx).name ': ' num2str([t u i]) ]);
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
