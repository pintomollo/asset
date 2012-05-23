function [chi2ple, psple] = model_identifiability(temp)

  if (nargin == 0)
    temp = 10;
  end

  opts = get_struct('modeling');
  opts = load_parameters(opts, 'goehring.txt');

  x0 = opts.init_params;
  x0 = repmat(x0, [opts.nparticles, 1]);

  flow = opts.advection_params;
  if (size(flow, 1) ~= size(x0, 1))
    [X, Y] = meshgrid([1:size(flow, 1)], 1+([0:size(x0, 1)-1]*(size(flow, 1)-1)/(size(x0, 1)-1)).');
    flow = bilinear_mex(flow, X, Y, [2 2]);
  end

  ml_params = [opts.diffusion_params; ...
                opts.reaction_params(1:end-2, :)];

  [orig, orig_t] = simulate_model(x0, [ml_params; opts.reaction_params(end-1:end, :)], opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow, opts.user_data);
  orig = orig((end/2)+1:end, :);

  penalty = -((3*max(orig(:)))^2)*opts.nparticles/2;

  rescaling = 10.^(floor(log10(ml_params)));
  rescaling(ml_params == 0) = 1;
  ml_params = ml_params ./ rescaling;

  uuid = now + cputime;
  RandStream.setDefaultStream(RandStream('mt19937ar','Seed',uuid));

  for i=1:numel(ml_params)
    [chi2ple{i}, psple{i}] = ple(ml_params(:), i, 20, @chi2score);
  end

  save(['ple-' num2str(uuid) '.mat'], 'chi2ple', 'psple');

  return;

  %function L = likelihood(params, stepping)
  function L = chi2score(params)

    %if (nargin == 1)
    %  stepping = 1;
    %end

    %[res, t] = simulate_model(x0, [params .* rescaling; opts.reaction_params(end-1:end, :)], opts.x_step, opts.tmax, opts.time_step / stepping, opts.output_rate, flow, opts.user_data);
    params = reshape(params, size(rescaling));
    [res, t] = simulate_model(x0, [params .* rescaling; opts.reaction_params(end-1:end, :)], opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow, opts.user_data);
    res = res((end/2)+1:end, :);
    if (length(t) ~= length(orig_t))
      res = interp1q(t.', res.', orig_t.').';
    end

    L = sum(((orig - res) / temp).^2);
    L(isnan(L)) = penalty;
    L = sum(L);
    %if (isnan(L))
    %  stepping = 2*stepping;
    %  if (stepping > 8)
    %    L = -Inf;
    %  else
    %    L = likelihood(params, stepping);
    %  end
    %end

    return;
  end
end
