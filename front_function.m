function [values, arctan] = front_function(params, pos)

  if (nargin == 1)
    npos = 16;
    pos = [-npos:npos];
  end

  if (numel(params) > 6)
    nframes = size(params, 1);
    values = NaN(nframes, length(pos));

    for i=1:nframes
      values(i, :) = front_function(params(i, :), pos);
    end

    return;
  end

  params = real(params);

  if (any(isinf(params) | isnan(params)))
    values = zeros(size(pos));

    return;
  end

  center = params(1);
  sigma = params(2);
  ampl = params(3);
  bkg = params(4);
  slope = params(5);
  atan_ampl = params(6);

  values = (atan((pos - center) * slope)/pi + 0.5) * atan_ampl + bkg; 
  if (nargout == 2)
    arctan = values;
    values = exp(-((pos - center).^2) / (2*(sigma^2)))*ampl;
  else
    values = values + exp(-((pos - center).^2) / (2*(sigma^2)))*ampl;
  end

  %figure;subplot(211);
  %plot((atan((pos - center) * slope)/pi + 0.5) * atan_ampl); 
  %subplot(212)
  %plot(exp(-((pos - center).^2) / (2*(sigma^2)))*ampl);

  return;
end
