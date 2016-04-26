function ml_values = extract_model_parameters(ml_values, convert_params)

  if (nargin < 2)
    convert_params = false;
  end

  opts = get_struct('modeling');
  opts = load_parameters(opts, 'goehring.txt');
  opts = load_parameters(opts, 'extended_model.txt');

  date_change = datenum(2013, 12, 01);
  params_struct = get_struct('model_parameter');

  for i = 1:size(ml_values, 1)
    fitting = ml_values{i, 1}{1};

    if (strncmp(fitting.fitting_type, 'sample', 6))
      norm_func = @(x)(x);
    else
      norm_func = @abs;
    end

    [fit_params, fit_energy, fit_temperatures, fit_viscosity] = model_description(fitting.parameter_set);

    nparams = length(fit_params);
    ngroups = max(length(fitting.temperature), 1);
    noffsets = ngroups*strncmp(fitting.aligning_type, 'fitting', 7);

    [temperatures, indx_fwd, indx_bwd] = unique(fitting.temperature);
    good_visc = (temperatures ~= opts.reaction_temperature);
    bad_flow = (temperatures == opts.flow_temperature);

    nvisc = sum(good_visc);
    ntemps = length(temperatures);
    nenergy = length(fit_energy);

    viscosities = ones(1, ntemps);

    for j = 1:size(ml_values{i, 2}, 1)
      value = ml_values{i, 2}{j, 2}(end);

      data = params_struct;
      effect_params = params_struct;

      ntotal = 0;

      if (~isfinite(value.score))
        value.score = Inf;
        warning('Empty optimization file !');
      else

        if (isempty(value.evolution))
          pts = [value.score value.params];
        else
          pts = value.evolution{1};
        end

        [npts, nvals] = size(pts);

        if (~isempty(fitting.fixed_parameter))
          if (size(pts, 2) ~= sum(~fitting.fixed_parameter)+1)
            value.score = Inf;
            warning('Warning: wrong number of parameters');
          else
            nfixed = length(fitting.fixed_parameter);

            if (nfixed == length(fitting.init_pos))
              tmp_pts = repmat(fitting.init_pos, npts, 1);
              tmp_pts(:, ~fitting.fixed_parameter) = pts(:,2:end);
            elseif (nfixed == length(fitting.init_pos) + fitting.scale_flow)
              tmp_pts = repmat([fitting.init_pos 1], npts, 1);
              tmp_pts(:, ~fitting.fixed_parameter) = pts(:,2:end);
            elseif (nfixed == length(fitting.init_pos) + noffsets)
              tmp_pts = NaN(npts, nfixed);
              if (nparams > 0)
                tmp_pts(:,1:nparams) = repmat(fitting.init_pos(:,1:nparams), npts, 1);
              end
              tmp_pts(:,nparams+noffsets+1:end) = repmat(fitting.init_pos(:,nparams+1:end), npts, 1);
              tmp_pts(:, ~fitting.fixed_parameter) = pts(:,2:end);
            elseif (nfixed == length(fitting.init_pos) + noffsets + fitting.scale_flow)
              tmp_pts = NaN(npts, nfixed);
              if (nparams > 0)
                tmp_pts(:,1:nparams) = repmat(fitting.init_pos(:,1:nparams), npts, 1);
              end
              tmp_pts(:,nparams+noffsets+1:end) = repmat([fitting.init_pos(:,nparams+1:end) opts.scale_params(1)], npts, 1);
              tmp_pts(:, ~fitting.fixed_parameter) = pts(:,2:end);
            elseif (nfixed == length(fitting.init_pos) + nparams + fitting.scale_flow + fitting.fit_flow)
              tmp_pts = NaN(npts, nfixed);
              tmp_mut = reshape(fitting.simulation_parameters, [], 2);
              tmp_pts(:,1:nparams) = repmat(tmp_mut(fit_params), npts, 1);
              tmp_pts(:,nparams+1:nparams+noffsets) = repmat(fitting.init_pos, npts, 1);
              if (fitting.fit_flow)
                tmp_pts(:,end-1) = 1;
              end
              if (fitting.scale_flow)
                tmp_pts(:,end) = 1;
              end
              tmp_pts(:, ~fitting.fixed_parameter) = pts(:,2:end);
            else

              tmp_pts = NaN(npts, nfixed);
              tmp_mut = reshape(fitting.simulation_parameters, [], 2);
              tmp_pts(:,1:nparams) = repmat(tmp_mut(fit_params), npts, 1);
              if (fitting.scale_flow)
                tmp_pts(:,end) = 1;
              end
              tmp_pts(:, ~fitting.fixed_parameter) = pts(:,2:end);
            end

            pts = [pts(:,1), tmp_pts];
            nvals = size(pts, 2);
          end
        end

        if (nvals > nparams)

          data.score = pts(:, 1);
          if (nparams > 0)
            data.rate = bsxfun(@times, norm_func(pts(:, 2:(nparams+1))), fitting.rescale_factor);
            ntotal = ntotal + nparams;
          end

          [good_params, params_indx] = ismember([4 5 12 13], fit_params);

          if (convert_params && fitting.fit_relative && all(good_params))
            data.rate(:, params_indx(1)) = data.rate(:, params_indx(1)) .* data.rate(:, params_indx(4));
            data.rate(:, params_indx(3)) = data.rate(:, params_indx(3)) .* data.rate(:, params_indx(2));
          end

          pts = pts(:, (nparams+2):end);
          data.offset = pts(:, 1:noffsets)*fitting.offset_scaling;
          ntotal = ntotal + noffsets;

          if (~isempty(pts))
            is_deprecated = (date_change - datenum(datevec(value.time)) > 0);

            if (fitting.scale_flow)
              data.flow_scaling = norm_func(pts(:, end));
              pts = pts(:, 1:end-1);
              ntotal = ntotal + 1;
            end

            if (fitting.fit_sigma)
              data.sigma = norm_func(pts(:, end));
              pts = pts(:, 1:end-1);
              ntotal = ntotal + 1;
            end

            if (fitting.fit_flow)
              data.flow = norm_func(pts(:, end));
              pts = pts(:, 1:end-1);
              ntotal = ntotal + 1;
            end

            if (fit_viscosity)
              data.viscosity = ones(npts, ntemps);
              if (is_deprecated)
                data.viscosity(:, good_visc(:, indx_bwd)) = norm_func(pts(:, end-nvisc+1:end));
                data.viscosity = data.viscosity(:, indx_fwd);
              else
                data.viscosity(:, good_visc) = norm_func(pts(:, end-nvisc+1:end));
              end
              data.viscosity = data.viscosity(:, good_visc);

              pts = pts(:, 1:end-nvisc);
              ntotal = ntotal + nvisc;
            end

            if (fit_temperatures && nenergy > 0)
              if (fitting.fit_model)
                data.energy = norm_func(pts(:,end-nenergy+1:end));
                ntotal = ntotal + nenergy;
              else
                data.energy = norm_func(pts(:,end-nenergy*nvisc+1:end));
                ntotal = ntotal + nenergy*nvisc;
              end
            end
          end
        end

        if (npts > 1)
          [junk, indx] = min(data.score);
          value.evolution = data;

          fields = fieldnames(data);

          for f = 1:length(fields)
            if (~isempty(data.(fields{f})))
              data.(fields{f}) = data.(fields{f})(indx, :);
            end
          end
        end

        effect_params.rate = ones(1, ntemps);
        effect_params.flow = ones(1, ntemps);
        effect_params.viscosity = ones(1, ntemps);

        if (fit_viscosity)
          effect_params.viscosity(good_visc) = data.viscosity;
        end

        if (fit_temperatures)
          if (fitting.fit_model)
            kB = 8.6173324e-5;
            C2K = 273.15;

            visc_E = 0;
            if (nenergy > 1)
              E = data.energy(1:end-1);
              E = E(:);
              flow_E = data.energy(end);
              if (nenergy == 3)
                visc_E = flow_E;
                flow_E = E(end);
                E = E(1);
              end
            elseif (nenergy > 0)
              E = data.energy;
              flow_E = data.energy;
            else
              E = 0.65;
              flow_E = E;
            end

            E = abs(E);
            flow_E = abs(flow_E);

            diff_ratio = ((temperatures+C2K) / (opts.reaction_temperature+C2K)) .* ...
                         exp(-(visc_E/kB)*((1./(temperatures+C2K)) - (1/(opts.reaction_temperature+C2K))));
            rate_ratio = exp(bsxfun(@times, -(E/kB), ((1./(temperatures+C2K)) - (1/(opts.reaction_temperature+C2K)))));
            flow_ratio = exp(-(flow_E/kB).*((1./(temperatures+C2K)) - (1/(opts.flow_temperature+C2K))));
            diff_ratio = effect_params.viscosity .* diff_ratio;

          else
            curr_ratios = ones(nenergy,ntemps);
            curr_ratios(:,good_visc) = reshape(data.energy, nenergy, nvisc);

            if (nenergy > 1)
              rate_ratio = curr_ratios(1:end-1,:);
              flow_ratio = curr_ratios(end, :);
            else
              rate_ratio = curr_ratios;
              flow_ratio = curr_ratios;
            end

            if (fitting.fit_flow)
              flow_ratio(~good_visc) = data.flow;
            end

            diff_ratio = effect_params.viscosity;
          end

          effect_params.rate = abs(rate_ratio);
          effect_params.flow = abs(flow_ratio);
          effect_params.viscosity = abs(diff_ratio);

          effect_params.temperature = temperatures;
        end
      end

      data.effective_value = effect_params;

      data.nparams = ntotal;
      data.temperature = temperatures;
      value.params = data;

      ml_values{i,2}{j,2} = value;
    end
  end


  return;
end
