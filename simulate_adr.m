function [all_x, all_m] = simulate_adr(x0, opts)

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

  %opts.nparticles = 1024;

  if (isempty(x0))
    if (~isempty(opts.init_func) & isa(opts.init_func, 'function_handle'))
      x0 = opts.init_func(opts);
    else
      x0 = opts.init_params;
      x0 = repmat(x0, [opts.nparticles, 1]);
    end
  end
  %x0 = bsxfun(@times, [1:opts.nparticles].'/opts.nparticles, x0);

  all_params = [opts.diffusion_params; ...
                opts.reaction_params];

  flow = opts.advection_params;
  
  %opts.tmax = 30;
  %opts.output_rate = 1;
  %flow = [1:6].' * [1:0.5:4];
  %flow = zeros(3, 10);

  if (size(flow, 1) ~= size(x0, 1))
    [X, Y] = meshgrid([1:size(flow, 1)], 1+([0:size(x0, 1)-1]*(size(flow, 1)-1)/(size(x0, 1)-1)).');
    flow = bilinear_mex(flow, X, Y, [2 2]);
  end

  [all_m, t] = simulate_model(x0, all_params, opts.x_step, opts.tmax, opts.time_step, opts.output_rate, flow, opts.user_data, opts.max_iter);
  all_x = all_m;

  return;

  order = 2;
  scheme = 2;
  h = opts.x_step;
  pos_index = [0:opts.nparticles-1] * h;

  % Compute the corresponding times
  dt = opts.time_step;
  ntimes = opts.tmax / dt;
  t = [0:ntimes]*dt;

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
  %flow = opts.advection_params;

  %if (size(flow, 1) ~= size(x0, 1))
  %  [X, Y] = meshgrid(1+([0:size(x0, 2)-1]*(size(flow, 2)-1)/(size(x0, 2)-1)), [1:size(flow, 1)].');
  %  flow = bilinear_mex(flow, X, Y, [2 2]);
  %end

  index_flow = [1:size(flow, 1)].';
  time_flow = ones(size(index_flow));
  prev = [];
  prev_deriv = [];

  for i = 2:ntimes

    curr_flow = bilinear_mex(flow, 1 + (time_flow * t(i-1) / t_flow), index_flow, bounds);
    curr_flow(isnan(curr_flow)) = 0;

    dadvec = -finite_difference(bsxfun(@times, curr_flow, x), h, 1, order, 'central', fid);
    ddiff = bsxfun(@times, Ds, finite_difference(x, h, 2, order, 'central', fid));
    dreac = nates(x, h, params, fid);

    deriv = (dadvec + ddiff + dreac);
    if (scheme == 1)
      x = x + dt*deriv;
    else
      tmp_x = x + dt*deriv;

      curr_flow = bilinear_mex(flow, index_flow, time_flow * t(i+1) / t_flow, bounds);
      curr_flow(isnan(curr_flow)) = 0;

      dadvec = -finite_difference(bsxfun(@times, curr_flow, tmp_x), h, 1, order, 'central', fid);
      ddiff = bsxfun(@times, Ds, finite_difference(tmp_x, h, 2, order, 'central'));
      dreac = nates(tmp_x, h, params);

      x = x + dt*(deriv + dadvec + ddiff + dreac) / 2;
    end

    if (mod(t(i-1), opts.output_rate) >= mod(t(i), opts.output_rate))
      all_x(:, count) = x(:);
      count = count + 1;
    end
  end

  return;
end
