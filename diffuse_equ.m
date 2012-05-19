function du = diffuse_equ(u, h, order, boundary)

  switch boundary
    case 'periodic'
      switch order
        case 2
          du = (u(:, [end 1:end-1]) - 2*u + u(:, [2:end 1])) / h^2;
      end

    case 'symmetric'
      switch order
        case 2
          du = (u(:, [1 1:end-1]) - 2*u + u(:, [2:end end])) / h^2;
      end
  end

  return;
end
