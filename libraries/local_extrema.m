function [xmax,imax,xmin,imin] = local_extrema(x, nneigh)
% LOCAL_EXTREMA Gets the local extrema points from a time series.
%
%   [XMAX,IMAX,XMIN,IMIN] = LOCAL_EXTREMA(X) returns the local minima
%   and maxima points of the vector X ignoring NaN's, where
%     XMAX - maxima points
%     IMAX - indexes of the XMAX
%     XMIN - minima points
%     IMIN - indexes of the XMIN
%
%   [...] = LOCAL_EXTREMA(X, NNEIGHBORS) returns only the maximal/minimal
%   extrama values in a [-NNEIGHBORS +NNEIGHBORS] local window.
%
%   Based on EXTRAM by
%   Lic. on Physics Carlos Adrián Vargas Aguilera
%   Physical Oceanography MS candidate
%   UNIVERSIDAD DE GUADALAJARA 
%   Mexico, 2004
%
% Wilson lab, University of Otago
% Simon Blanchoud
% 25.06.2015

  % No local average
  if (nargin < 2)
    nneigh = 1;
  end

  xmax = [];
  imax = [];
  xmin = [];
  imin = [];

  % Vector input?
  Nt = numel(x);
  if Nt ~= length(x)
   error('Entry must be a vector.')
  end

  % NaN's:
  inan = find(isnan(x));
  indx = 1:Nt;
  if ~isempty(inan)
   indx(inan) = [];
   x(inan) = [];
   Nt = length(x);
  end

  % Difference between subsequent elements:
  x = x(:);
  dx = differentiator(x, 'super');

  % Local maxima and minima for windowing
  nneigh = ceil(nneigh);
  if (nneigh > 1)
    window = 2*nneigh + 1;
    indxs = [1:nneigh nneigh+2:window].';
    ssize = [window, 1];
    lmax = colfilt(x, ssize, 'sliding', @(y)(max(y(indxs,:), [], 1)));

    if (size(x, 1) == 1)
      lmax = lmax.';
    end

    % Only if needed
    if (nargout > 2)
      lmin = colfilt(x, ssize, 'sliding', @(y)(min(y(indxs,:), [], 1)));

      if (size(x, 1) == 1)
        lmin = lmin.';
      end
    end
  else
    lmax = x-1;
    lmin = x+1;
  end

  % Is an horizontal line?
  if ~any(dx)
   return
  end

  % Flat peaks? Put the middle element:
  a = find(dx~=0);              % Indexes where x changes
  lm = find(diff(a)~=1) + 1;    % Indexes where a do not changes
  d = a(lm) - a(lm-1);          % Number of elements in the flat peak
  a(lm) = a(lm) - floor(d/2);   % Save middle elements
  a(end+1) = Nt;

  % Peaks?
  xa  = x(a);             % Serie without flat peaks
  b = (diff(xa) > 0);     % 1  =>  positive slopes (minima begin)  
                          % 0  =>  negative slopes (maxima begin)
  xb  = diff(b);          % -1 =>  maxima indexes (but one) 
                          % +1 =>  minima indexes (but one)

  a = a(1:end-1);

  % Compensate for the missing datapoints after diff
  if (size(x,1)==1)
    xb = [0 xb];
  else
    xb = [0; xb];
  end

  imax = find((xb == -1) & (x(a) > lmax(a))); % maxima indexes
  imax = a(imax);

  % Only if needed
  if (nargout > 2)
    imin = find((xb == +1) & (x(a) < lmin(a))); % minima indexes
    imin = a(imin);
  else
    imin = [];
  end

  nmaxi = length(imax);
  nmini = length(imin);

  % Maximum or minumim on a flat peak at the ends?
  if (nmaxi==0) && (nmini==0)
   if x(1) > x(Nt)
    xmax = x(1);
    imax = indx(1);
    xmin = x(Nt);
    imin = indx(Nt);
   elseif x(1) < x(Nt)
    xmax = x(Nt);
    imax = indx(Nt);
    xmin = x(1);
    imin = indx(1);
   end
   return
  end
  xmax = x(imax);
  xmin = x(imin);

  % NaN's:
  if (nargout > 1)
    if ~isempty(inan)
     imax = indx(imax);
     imin = indx(imin);
    end

    % Same size as x:
    imax = reshape(imax,size(xmax));
    imin = reshape(imin,size(xmin));
  end

  return;
end
