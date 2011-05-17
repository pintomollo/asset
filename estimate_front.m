function params = estimate_front(signal, pos)

  wsize = floor(length(pos) / 3);
  
  inside = signal(1:wsize);
  outside = signal(end-wsize+1:end);
  cortex = signal(wsize+1:end-wsize);

  gap = [1:length(cortex)] / (length(cortex) + 1);
  arctan = (outside(1) - inside(end)) * gap + inside(end);
  arctan = [inside arctan outside];

  gauss = signal - arctan;
  gauss(gauss < 0) = 0;

  if (sum(gauss > 0) < 3)
    %arctan = (outside(1) - inside(end)) * gap + inside(end)*0.75;
    %arctan = [inside arctan outside];
    tmp = signal - arctan;
    tmp = tmp - min(tmp);
    gauss(wsize+1:end-wsize) = tmp(wsize+1:end-wsize);
  end

  ampl = max(gauss);
  gauss = gauss / sum(gauss);
  center = sum(gauss .* pos);
  sigma = sqrt(sum(gauss .* ((pos - center).^2)));
  values = exp(-((pos - center).^2) / (2*(sigma^2)))*ampl;

  arctan = signal - values;
  bkg = min(arctan);
  atan_ampl = max(arctan) - bkg;
  slope = diff(arctan) / atan_ampl;
  slope_pos = pos(1:end-1) - center;
  slope = ((1 + sqrt(1 - 4*(slope.^2).*(slope_pos.^2))) ./ (4*slope.*(slope_pos.^2)));
  slope = median(abs(slope) .* sign(real(slope)));

  %values = ((atan((pos - center) * slope)/pi + 0.5) * atan_ampl); 
  %values = values + exp(-((pos - center).^2) / (2*(sigma^2)))*ampl + bkg;

  %figure;plot(pos, signal);
  %hold on;
  %plot(pos, arctan, 'r')
  %plot(pos, signal - arctan, 'g')
  %plot(pos, values, 'k')

  params = [center, sigma, ampl, bkg, slope, atan_ampl];

  return;
end
