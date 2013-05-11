function [init_val, correct] = init_goehring(opts, is_relative)

  if (nargin==1)
    is_relative = false;
  end

  p = opts.reaction_params;
  correct = true;

  L = p(7);
  if (is_relative)
    p(5,:) = p(5,:) ./ (L^2);
    p(3,:) = p(3,:) .* p(4,:) .* L.^(2*p(4,:));
  end
  %figure;hold on;

  %[pts, lines] = study_isoclines(@maintenance, [0 10], opts.reaction_params);
  %figure;hold on;
  %plot(lines(:,1,1), lines(:,2,1), 'b');
  %plot(lines(:,1,2), lines(:,2,2), 'r');
  %scatter(pts(:,1), pts(:, 2), 'k');

  pos = p(7)*[0:(opts.nparticles-1)]/(opts.nparticles-1);
  full_pos = bsxfun(@minus, [pos; pos].', [0.45 0.55]*p(7));

  if (opts.init_params)

    pts = find_fixed_points(@uniform, [0 6], p);

    if (numel(pts) < 2)
      init_val = repmat([1.51 0.1], opts.nparticles, 1); % Values produced by the original parameters
      correct = false;
    else
      init_val = repmat(pts(1,:), opts.nparticles, 1);
    end

    %plot(pos, init_val(:,1), 'c');
    %plot(pos, init_val(:,2), 'm');

  else

    pts = find_fixed_points(@maintenance, [0 6], p);
    vals = max(pts);

    if (numel(vals) ~= 2)
      vals = [1.8 4]; % Values produced by the original parameters
      correct = false;
    end
    init_val = bsxfun(@times, vals, (erf(bsxfun(@times, [-2/9.8 2/7.7], full_pos)) + 1) / 2);
    %vals * (erf(a*(L-u)) + 1)/2;

    %plot(pos, init_val(:,1), 'b');
    %plot(pos, init_val(:,2), 'r');

  end

  if (is_relative)
    init_val / (L^2);
  end

  return;
end

function isos = uniform(pts, params)
% The function that defines the isoclines

  A = pts(:, 1);
  P = pts(:, 2);

  % The actual isoclines
  p = params(:,1);
  aiso = p(1)*p(5) ./ (p(1)*p(6) + p(2) + p(3)*P.^p(4));
  p = params(:,2);
  piso = p(1)*p(5) ./ (p(1)*p(6) + p(2) + p(3)*A.^p(4));

  % Concatenate them properly for study_isoclines
  isos = cat(3, [aiso P], [A piso]);

  return;
end

function isos = maintenance(pts, params)
% The function that defines the isoclines

  A = pts(:, 1);
  P = pts(:, 2);

  % The actual isoclines
  p = params(:,1);
  aiso = p(1)*p(5) ./ (p(1)*p(6)*(int_erf(2/9.8, p(7))+p(7))/(2*p(7)) + p(2) + p(3)*P.^p(4));
  p = params(:,2);
  piso = p(1)*p(5) ./ (p(1)*p(6)*(int_erf(2/7.7, p(7))+p(7))/(2*p(7)) + p(2) + p(3)*A.^p(4));

  % Concatenate them properly for study_isoclines
  isos = cat(3, [aiso P], [A piso]);

  return;
end

function int = int_erf(a, L)

  u = 0.55*L;
  int = (1/(sqrt(pi)*a))*(exp(-(a^2)*((L-u)^2)) - exp(-(a^2)*(u^2))) + (L-u)*erf(a*(L-u)) - u*erf(a*u);

  return;
end
