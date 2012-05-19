function du = finite_difference(u, dx, deriv, accuracy, method)
  
  switch method
    case 'forward'
      switch accuracy
        case 1
        case 2
      end
    case 'central'
      switch accuracy
        case 1
          if (deriv == 1)
            du = (-u(:, [1 1:end-1]) + u(:, [2:end end])) / (2*dx);
          else
            du = (u(:, [1 1:end-1]) - 2*u + u(:, [2:end end])) / dx^2;
          end
        case 2
          if (deriv == 1)
            du = (u(:, [2 1 1:end-2]) - 8*u(:, [1 1:end-1]) + 8*u(:, [2:end end]) - u(:, [3:end end end-1])) / (12*dx);
          else
            du = (-u(:, [2 1 1:end-2]) + 16*u(:, [1 1:end-1]) - 30*u + 16*u(:, [2:end end]) - u(:, [3:end end end-1])) / (12*(dx^2));
          end
      end
  end

  return;
end
