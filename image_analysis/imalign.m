function new_img = imalign(img, new_size, center, orientation)

  if (~isfloat(img))
    img = double(img);
  end

  if (size(new_size, 1) > 1)
    new_size = new_size.';
  end

  new_img = NaN(new_size); 
  new_center = new_size / 2;

  row = [1:new_size(2)].';
  row = row - new_center(2);
  col = ones(size(row));

  center = repmat(center(:), 1, length(row));

  rot_matrix = [cos(orientation) -sin(orientation); sin(orientation) cos(orientation)];

  for i=1:new_size(1)
    col_indx = (col * i) - new_center(1);
    pts = [row, col_indx] * rot_matrix; 
    pts = pts.' + center;

    new_img(i,:) = bilinear_mex(img, pts);
  end

  return;
end
