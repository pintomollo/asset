function cMap = purples(N)
% PURPLES creates a red and blue colormap as defined in [1].
%
%   PURPLES(M) returns an M-by-3 matrix containing a purple
%   sequential color palette. M is the number of different colors in the
%   colormap. If M is empty, a default value of 11 will be used.
%
% References:
%   [1] Brewer, Cynthia A., 2014. http://www.ColorBrewer.org, accessed date: 07.07.14
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 07.07.2014

  % As in other colormaps, allows to color existing axes
  if (nargin < 1 || isempty(N) || N == 0)
    N = 64;
  end

  % Get the colormap from brewer directly
  cMap = brewermap(N, 'Purples');

  return;
end
