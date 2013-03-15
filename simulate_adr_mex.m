function [all_x, t] = simulate_adr_mex(x0, opts)

  %%%% Adapt volume to surface ratio to real embryos

  if (nargin == 1)
    if (isstruct(x0))
      opts = x0;
      x0 = [];
    else
      opts = get_struct('modeling');
      opts = load_parameters(opts, 'goehring.txt');
    end
  elseif (nargin == 0)
    opts = get_struct('modeling');
    opts = load_parameters(opts, 'goehring.txt');
    x0 = [];
  end

  if (isempty(x0))
    if (~isempty(opts.init_func) & isa(opts.init_func, 'function_handle'))
      x0 = opts.init_func(opts);
    else
      x0 = opts.init_params;
      x0 = repmat(x0, [opts.nparticles, 1]);
    end
  end

  all_params = [opts.diffusion_params; ...
                opts.reaction_params];

  flow = opts.advection_params;
  if (size(flow, 1) ~= size(x0, 1))
    [X, Y] = meshgrid([1:size(flow, 2)], 1+([0:size(x0, 1)-1]*(size(flow, 1)-1)/(size(x0, 1)-1)).');
    flow = bilinear_mex(flow, X, Y, [2 2]);
  end

  [all_m, t] = simulate_model(x0, all_params, opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow, opts.user_data, opts.max_iter);
  all_x = all_m;

  return;
end
