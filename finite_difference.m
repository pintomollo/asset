function du = finite_difference(u, dx, deriv, accuracy, method, fid)
  
  switch method
    case 'forward'
      %switch accuracy
      %  case 1
          du = (u - u([1 1:end-1], :)) / dx;
          dun = (-u + u([2:end end], :)) / dx;
          du(u < 0) = dun(u < 0);
      %  case 2
      %end
    case 'central'
      switch accuracy
        case 1
          if (deriv == 1)
            du = (-u([1 1:end-1], :) + u([2:end end], :)) / (2*dx);
          else
            du = (u([1 1:end-1], :) - 2*u + u([2:end end], :)) / dx^2;
          end
        case 2
          if (deriv == 1)
            du = (u([2 1 1:end-2], :) - 8*u([1 1:end-1], :) + 8*u([2:end end], :) - u( [3:end end end-1], :)) / (12*dx);
%            if (nargin == 6)
%              fprintf(fid, '%f (%d) %f %f %f %f %f\n', [du(1,:); [1:size(u, 2)]; u(1, [2 1 1:end-2]); u(1, [1 1:end-1]); u(1,:); u(1, [2:end end]); u(1, [3:end end end-1])]);
%              fprintf(fid, '...........\n');
%            end

          else
            du = (-u([2 1 1:end-2], :) + 16*u([1 1:end-1], :) - 30*u + 16*u([2:end end],:) - u([3:end end end-1],:)) / (12*(dx^2));
%            if (nargin == 6)
%              fprintf(fid, '%f (%d) %f %f %f %f %f\n', [du(1,:); [1:size(u, 2)]; u(1, [2 1 1:end-2]); u(1, [1 1:end-1]); u(1,:); u(1, [2:end end]); u(1, [3:end end end-1])]);
%              fprintf(fid, '...........\n');
%            end
          end
      end
  end
  
  return;
end
