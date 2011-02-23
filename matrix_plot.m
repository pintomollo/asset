function hhh = matrix_plot(mat)

  if(~isnumeric(mat))
    mat = double(mat);
  end

  h = gca;
  hh = imagesc(mat);
  set(hh, 'Parent', h);
  set(h, 'YDir', 'normal');

  if(nargout>0)
    hhh=hh;
  end

  return;
end
