function perif = ellipse_circum(axes_length, aim)

  if (nargin == 2)
    perif = fzero(@(u)(ellipse_circum(axes_length*u)-aim), [0.25 4]);
  else
    a = axes_length(1,:);
    b = axes_length(2,:);

    ratio = 3*(((a - b) ./ (a + b)).^2);
    perif = pi*(a+b).*(1+(ratio./(10+sqrt(4-ratio))));

    perif(~isfinite(perif)) = 0;
  end

  return;
end
