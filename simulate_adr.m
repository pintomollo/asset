function [all_x, t] = simulate_adr(x0, opts)

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
    [X, Y] = meshgrid([1:size(flow, 1)], 1+([0:size(x0, 1)-1]*(size(flow, 1)-1)/(size(x0, 1)-1)).');
    flow = bilinear_mex(flow, X, Y, [2 2]);
  end

  order = 2;
  scheme = 1;
  h = opts.x_step;

  max_iter = opts.max_iter;

  % Compute the corresponding times
  dt = opts.time_step;
  ntimes = opts.tmax / dt;

  % Initialize the output
  all_x = zeros(numel(x0), ceil(opts.tmax/opts.output_rate));
  all_x(:, 1) = x0(:);
  x = x0;

  count = 2;

  bounds = [2 3];
  t_flow = opts.user_data;
  fid = [];

  Ds = opts.diffusion_params;
  params = opts.reaction_params;

  [npts, nflow] = size(flow);
  prev = [];
  prev_deriv = [];

  total_time = 0;
  time_count = 0;

  [t, all_x] = ode45(@deriv_func, [0 opts.tmax], x);

  return;

  function dy = deriv_func(t, y)
    
    y = reshape(y, npts, []);

    time_flow = t/t_flow;
    lindx = floor(time_flow) + 1;

    if (lindx >= nflow)
      curr_flow = zeros(npts, 1);
    else
      curr_flow = flow(:, lindx)*(1-(time_flow - lindx)) + flow(:, lindx+1)*(time_flow - lindx);
    end

    dadvec = -finite_difference(bsxfun(@times, curr_flow, y), h, 1, order, 'forward', fid);
    ddiff = bsxfun(@times, Ds, finite_difference(y, h, 2, 2, 'central', fid));
    dreac = nates(y, h, params, fid);

    dy = (dadvec + ddiff + dreac);
    dy = dy(:);

    return;
  end
end
