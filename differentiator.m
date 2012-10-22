function dy = differentiator(varargin)
% DIFFERENTIATOR computes the numerical derivative of a provided function using various
% methods.
%
%   DY = DIFFERENTIATOR(X, Y, DIM, N, METHOD) computes the first derivative DY for the
%   function Y = f(X). For multi-dimensional Y, it computes it along the dimension DIM.
%   The numerical derivative is computed using METHOD with a neighborhood of N. Available
%   methods are:
%     - 'cfd'     : Central Finite Difference (Polynomial approximation of the signal)
%     - 'lanczos' : Low-noise Lanczos (a special Savitzky-Golay digital differentiator,
%                   CFD with least-square approximation to reduce Gaussian noise)
%     - 'noise'   : Smooth noise-robust (smooth Lanczos used when the noise is
%                   restricted to a certain range of frequencies)
%
%   The implemented neighborhoods N depend on each method but typically range from 5 to
%   11 (all these methods are central and thus have an odd neighborhood).
%
%   DY = DIFFERENTIATOR(..., IS_SUPER) uses the 'super' differentiator instead of the
%   standard one, thus using polynomes of degree 4 instead of 2 for the noise supression.
%   This option is not available for the CFD method.
%
%   DY = DIFFERENTIATOR(..., BOUNDARY) specifies the type of boundary condition (necessary
%   to obtain DY with as many elements as Y. Available types of boundary condition are:
%     - 'circular'  : consider Y as being circular (Y(1) is next to Y(end))
%     - 'replicate' : repeats the border elements
%     - 'symmetric' : uses a mirror reflection of Y
%
%   DY = DIFFERENTIATOR(Y) computes the derivative assuming:
%     - X        = uniformly-spaced points
%     - DIM      = first non-singleton dimension
%     - N        = 7
%     - METHOD   = 'lanczos'
%     - IS_SUPER = false
%     - BOUNDARY = 'symmetric'
%
% Both the methods and the coefficients have been developed by Pavel Holoborodko:
% http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 16.12.2011

  [x, y, dim, nneigh, method, super, boundary] = parse_input(varargin{:});

  if (isempty(y))
    dy = y;

    return;
  end

  coefs = get_coefs(method, nneigh, super);
  nrepeat = length(coefs);
  center = nrepeat + 1;

  sizey = size(y);
  perm_dim = [1:length(sizey)];
  perm_dim(dim) = 1;
  perm_dim(1) = dim;

  y = permute(y, perm_dim);
  x = permute(x, perm_dim);

  y = reshape(y, sizey(dim), []);
  x = reshape(x, sizey(dim), []);

  dy = NaN(size(y));

  valids = all(isfinite(y), 2) & all(isfinite(x), 2);
  y = y(valids, :);
  x = x(valids, :);

  index = [1:size(y,1)] + center - 1;

  y = padarray(y, nrepeat, boundary, 'both');
  x = padarray(x, nrepeat, 'symmetric', 'both');
  x(1:nrepeat, :) = bsxfun(@minus, 2*x(nrepeat + 1, :), x(1:nrepeat, :)) - 1;
  x(end-nrepeat+1:end, :) = bsxfun(@minus, 2*x(end-nrepeat, :), x(end-nrepeat+1:end, :)) + 1;

  if (~isempty(y))
    accum = zeros(size(y) - [2*nrepeat 0]);

    for k=1:nrepeat
      accum = accum + 2*k*coefs(k)* bsxfun(@rdivide, y(index+k, :) - y(index-k, :), x(index+k, :) - x(index - k, :));
    end
  else
    accum = [];
  end

  dy(valids, :) = accum;

  dy = reshape(dy, sizey(perm_dim));
  dy = ipermute(dy, perm_dim);

  return;
end

function coefs = get_coefs(method, nneigh, super)

  switch method
    case 'cfd'
      switch nneigh
        case 3
          coefs = [0.5];
        case 5
          coefs = [8 -1]/12;
        case 7
          coefs = [45 -9 1]/60;
        case 9
          coefs = [672 -168 32 -3]/840;
        otherwise
          coefs = [672 -168 32 -3]/840;
          warning(['The differentiator ' method ' with N=' num2str(nneigh) ' is not implemented, using N=9 instead']);
      end
    case 'lanczos'
      switch nneigh
        case 5
          coefs = [1 2]/10;
        case 7
          if (super)
            coefs = [58 67 -22]/252;
          else
            coefs = [1 2 3]/28;
          end
        case 9
          if (super)
            coefs = [126 193 142 -86]/1188;
          else
            coefs = [1 2 3 4]/60;
          end
        case 11
          if (super)
            coefs = [296 503 532 294 -300]/5148;
          else
            coefs = [1 2 3 4 5]/110;
          end
        otherwise
          if (super)
            coefs = [296 503 532 294 -300]/5148;
            warning(['The differentiator super-' method ' with N=' num2str(nneigh) ' is not implemented, using N=11 instead']);
          else
            coefs = 12*[1:((nneigh-1)/2)]/(nneigh^3 - nneigh);
          end
      end
    case 'noise'
      switch nneigh
        case 5
          coefs = [2 1]/8;
        case 7
          if (super)
            coefs = [39 12 -5]/96;
          else
            coefs = [5 4 1]/32;
          end
        case 9
          if (super)
            coefs = [27 16 -1 -2]/60;
          else
            coefs = [14 14 6 1]/128;
          end
        case 11
          if (super)
            coefs = [322 256 39 -32 -11]/1536;
          else
            coefs = [42 48 27 8 1]/512;
          end
        otherwise
          if (super)
            coefs = [322 256 39 -32 -11]/1536;
          else
            coefs = [42 48 27 8 1]/512;
          end
          warning(['The differentiator ' method ' with N=' num2str(nneigh) ' is not implemented, using N=11 instead']);
      end
  end
  
  return;
end

function [x, y, dim, nneigh, method, super, boundary] = parse_input(varargin)

  x = [];
  y = [];
  dim = 0;
  nneigh = 7;
  method = 'lanczos';
  boundary = 'symmetric';
  super = false;

  for i = 1:length(varargin)
    var_type = get_type(varargin{i});
    switch var_type
      case 'none'
        dim = 1;
      case 'num'
        if (numel(varargin{i}) == 1)
          if (dim == 0)
            dim = varargin{i};
          else
            nneigh = varargin{i};
          end
        elseif (isempty(x))
          x = varargin{i};
        elseif (isempty(y))
          y = varargin{i};
        end
      case 'bool'
        super = varargin{i};
      case 'char'
        switch varargin{i}
          case {'symmetric', 'replicate', 'circular'}
            boundary = varargin{i};
          case 'super'
            super = true;
          otherwise
            method = varargin{i};
        end
    end
  end

  if (isempty(y))
    y = x;
    x = [];
  end

  npts = size(y);
  good_dims = (npts > 1);
  first_dim = find(good_dims, 1, 'first');

  if (~isempty(first_dim))
    if (dim == 0 | ~good_dims(dim))
      dim = first_dim;
    end

    if (isempty(x))
      tmp_size = ones(size(npts));
      tmp_size(dim) = npts(dim);
      x = reshape([1:npts(dim)], tmp_size);
    end
  end

  if (size(x, dim) ~= npts(dim))
    error(['The dimensions of x (' num2str(size(x,dim)) ') and y (' num2str(npts(dim)) ') do not correspond.']);
  end

  nneigh = nneigh + 1 - mod(nneigh, 2);

  switch method
    case 'cfd'
      if (nneigh < 3)
        nneigh = 3;
      end
    case {'lanczos', 'noise'}
      if (super & nneigh < 7)
        nneigh = 7;
      elseif (~super & nneigh < 5)
        nneigh = 5;
      end
  end

  return;
end
