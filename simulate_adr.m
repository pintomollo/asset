function [all_x] = simulate_adr(x0, opts, order, scheme)

  %keyboard
  if (nargin < 4)
    order = 2;
    scheme = 1;
  end

  h = diff(opts.boundaries) / opts.nparticles;
  all_x = [];

  t_flow = opts.user_data.time_steps;
  dt = opts.time_step;
  tmax = opts.tmax;
  output_rate = opts.output_rate;

  all_params = [opts.diffusion_params; ...
                [[opts.reaction_params(1:3) 1 opts.reaction_params([4 5 11 12])].' ...
                [opts.reaction_params(6:9) 1 opts.reaction_params([10 11 12])].']];

  flow = opts.advection_params;

  if (size(flow, 2) ~= size(x0, 2))
    [X, Y] = meshgrid([1:size(x0, 2)]*size(flow, 2)/size(x0, 2), [1:size(flow, 1)].');
    flow = bilinear_mex(flow, X, Y, [2 2]);
  end

  all_x = simulate_model(x0.', all_params, h, tmax, dt, output_rate, flow.', t_flow);

  return;

  pos_index = [0:opts.nparticles-1] * h;

  % Compute the corresponding times
  dt = opts.time_step;
  ntimes = opts.tmax / dt;
  t = [0:ntimes]*dt;

  % Initialize the output
  all_x = zeros([size(x0) ceil(opts.tmax/opts.output_rate)]);
  all_x(:, :, 1) = x0;
  x = x0;

  count = 2;

  bounds = opts.user_data.boundary_type;
  t_flow = opts.user_data.time_steps;

  Ds = opts.diffusion_params.';
  params = opts.reaction_params;
  flow = opts.advection_params;

  if (size(flow, 2) ~= size(x0, 2))
    [X, Y] = meshgrid([1:size(x0, 2)]*size(flow, 2)/size(x0, 2), [1:size(flow, 1)].');
    flow = bilinear_mex(flow, X, Y, [2 2]);
  end

  index_flow = [1:size(flow, 2)];
  time_flow = ones(size(index_flow));
  prev = [];
  prev_deriv = [];

  for i = 2:ntimes

    curr_flow = bilinear_mex(flow, index_flow, time_flow * t(i) / t_flow, bounds);
    curr_flow(isnan(curr_flow)) = 0;

    dadvec = -finite_difference(bsxfun(@times, curr_flow, x), h, 1, order, 'central')
    ddiff = bsxfun(@times, Ds, finite_difference(x, h, 2, order, 'central'))
    dreac = nates(x, h, params)

    %x = x + dt*(dadvec + ddiff + dreac);
    deriv = (dadvec + ddiff + dreac);
    if (scheme == 1)
      x = x + dt*deriv
    else
      tmp_x = x + dt*deriv;

      curr_flow = bilinear_mex(flow, index_flow, time_flow * t(i+1) / t_flow, bounds);
      curr_flow(isnan(curr_flow)) = 0;

      dadvec = -finite_difference(bsxfun(@times, curr_flow, tmp_x), h, 1, order, 'central');
      ddiff = bsxfun(@times, Ds, finite_difference(tmp_x, h, 2, order, 'central'));
      dreac = nates(tmp_x, h, params);

      x = x + dt*(deriv + dadvec + ddiff + dreac) / 2;
    end

    if (mod(t(i-1), opts.output_rate) >= mod(t(i), opts.output_rate))
      all_x(:,:,count) = x;
      count = count + 1;

      %plot(pos_index, dadvec(1,:), 'b');
      %title(num2str(t(i)))
      %drawnow
    end
    
    break;
  end

  return;
end
