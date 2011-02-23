function perif = ellipse_circum(axes_length)

  a=axes_length(1,:);
  b=axes_length(2,:);

  perif = pi*(a+b).*(1+(3*((a-b)./(a+b)).^2)./(10+sqrt(4-3*((a-b)./(a+b)).^2)));

  if (isnan(perif) | isinf(perif))
    perif = 0;
  end

  return;
end
