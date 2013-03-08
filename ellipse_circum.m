function perif = ellipse_circum(axes_length)

  a = axes_length(1,:);
  b = axes_length(2,:);

  ratio = 3*(((a - b) ./ (a + b)).^2);
  perif = pi*(a+b).*(1+(ratio./(10+sqrt(4-ratio))));

  perif(~isfinite(perif)) = 0;

  return;
end
