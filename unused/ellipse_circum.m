function perif = ellipse_circum(axes_length, aim, x_only)

  if (nargin < 3)
    x_only = false;
  end

  permute = false;
  if (nargin > 1 && size(aim, 1) > 1)
    aim = aim.';
    axes_length = axes_length.';
    permute = true;
  end

  if (nargin ~= 1)
    perif = NaN(1, length(aim));

    for i=1:length(aim)
      if (x_only)
        tmp_axes = axes_length(:, i);
        perif(i) = fzero(@(u)(ellipse_circum([u; tmp_axes(2:3)])-aim(i)), [0.25 4]*tmp_axes(1));
      else
        perif(i) = fzero(@(u)(ellipse_circum(axes_length(:,i)*u)-aim(i)), [0.25 4]);
      end
    end
  else
    a = axes_length(1,:);
    b = axes_length(2,:);

    ratio = 3*(((a - b) ./ (a + b)).^2);
    perif = pi*(a+b).*(1+(ratio./(10+sqrt(4-ratio))));

    perif(~isfinite(perif)) = 0;
  end

  if (permute)
    perif = perif.';
  end

  return;
end
