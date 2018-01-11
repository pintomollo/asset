function [ratio, surface, volume] = surface2volume(axes_length)

  axes_length = sort(axes_length, 'descend');

  a = axes_length(1,:);
  b = axes_length(2,:);
  c = axes_length(3,:);

  npts = numel(a);
  if (npts > 1)
    ratio = NaN(1, npts);
    surface = ratio;
    volume = ratio;

    for i=1:npts
      [ratio(i), surface(i), volume(i)] = surface2volume(axes_length(:,i));
    end

    return;
  end

  phi = acos(c./a);
  K = (a.^2 .* (b.^2 - c.^2)) ./ (b.^2 .* (a.^2 - c.^2));
  bads = ~isfinite(K);
  K(bads) = 0;
  [F, E] = elliptic12(phi, K);

  surface = 2.*pi.*c.^2 + ((2.*pi.*a.*b)./(sin(phi))).*(E.*sin(phi).^2 + F.*cos(phi).^2);
  volume = 4/3 .* pi .* a .* b .* c;

  ratio = surface ./ volume;
  ratio(bads) = 3 ./ a(bads);

  return;
end
