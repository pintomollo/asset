function perif = ellipse_circum(axes_length, aim, x_only)

  if (nargin < 3)
    x_only = false;
  end

  if (nargin ~= 1)
    if (x_only)
      axes_length = axes_length(:);
      perif = fzero(@(u)(ellipse_circum([u; axes_length(2:3)])-aim), [0.25 4]*axes_length(1));
    else
      perif = fzero(@(u)(ellipse_circum(axes_length*u)-aim), [0.25 4]);
    end
  else
    a = axes_length(1,:);
    b = axes_length(2,:);

    ratio = 3*(((a - b) ./ (a + b)).^2);
    perif = pi*(a+b).*(1+(ratio./(10+sqrt(4-ratio))));

    perif(~isfinite(perif)) = 0;
  end

  return;
end
