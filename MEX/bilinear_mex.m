% BILINEAR_MEX compute a bilinear interpolation of an image at the provided
% (sub)pixel coordinates.
%
%   PIXS = BILINEAR_MEX(IMG, XCOORD, YCOORD) returns the bilinear interpolation of the
%   pixel values PIXS from IMG at the provided coordinates (XCOORD, YCOORD) tuples.
%   PIXS has the same dimensions as XCOORD. Note that coordinates should be provided
%   as carthesian coordinates, not as matrix indexes !
%
%   PIXS = BILINEAR_MEX(IMG, COORDS) where COORDS is a Nx2 matrix with each row
%   corresponding to a [X_coord, X_coord] tuple.
%
%   PIXS = BILINEAR_MEX(..., BOUNDARY) defines in addition the type of boundary
%   conditions to be used. Available behaviors are: 0 (NaN outside), 1 (circular),
%   2 (replicate), 3 (symmteric). See imfilter for more details on such behaviors.
%   If BOUNDARY has two elements, X and Y behaviors can be defined separately.
%   By default, both coordinates are set to 0.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 07.07.2014
