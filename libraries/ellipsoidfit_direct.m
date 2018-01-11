function [center, radii, R] = ellipsoidfit_direct(x,y,z)
% Direct least squares fitting of ellipsoids under the constraint 4J - I^2 > 0.
% The constraint confines the class of ellipsoids to fit to those whose smallest radius
% is at least half of the largest radius.
%
% Input arguments:
% x,y,z;
%    x, y and z coodinates of 3D points
%
% Output arguments:
% center:
%    ellispoid center coordinates [cx; cy; cz]
% ax:
%    ellipsoid semi-axes (radii) [a; b; c]
% R:
%    ellipsoid rotation (radii directions as rows of the 3x3 matrix)
%
% References:
% Qingde Li and John G. Griffiths, "Least Squares Ellipsoid Specific Fitting",
%    Proceedings of the Geometric Modeling and Processing, 2004.

% Copyright 2011 Levente Hunyadi

narginchk(3,3);
validateattributes(x, {'numeric'}, {'real','nonempty','vector'});
validateattributes(y, {'numeric'}, {'real','nonempty','vector'});
validateattributes(z, {'numeric'}, {'real','nonempty','vector'});
x = x(:);
y = y(:);
z = z(:);

% build design matrix
D = [ x.^2, y.^2, z.^2, 2*y.*z, 2*x.*z, 2*x.*y, 2*x, 2*y, 2*z, ones(numel(x),1) ];

% build scatter matrix
S = D'*D;

% build 10x10 constraint matrix
k = 4;  % to ensure that the parameter vector always defines an ellipse
C1 = [ 0 k k ; k 0 k ; k k 0 ] / 2 - 1;
C2 = -k * eye(3,3);
C3 = zeros(4,4);

C = blkdiag(C1,C2,C3);

% solve eigensystem
[gevec, geval] = eig(S,C);
geval = diag(geval);

% extract eigenvector corresponding to the unique positive eigenvalue
flt = geval > 0 & ~isinf(geval);
switch nnz(flt)
    case 0
        % degenerate case; single positive eigenvalue becomes near-zero negative eigenvalue
        % due to round-off error
        [~,ix] = min(abs(geval));
        v = gevec(:,ix);
    case 1
        % regular case
        v = gevec(:,flt);
    otherwise
        % degenerate case; several positive eigenvalues appear
        [~,ix] = min(abs(geval));
        v = gevec(:,ix);
end

p = zeros(size(v));
p(1:3) = v(1:3);
p(4:6) = 2*v(6:-1:4);  % exchange order of y*z, x*z, x*y to x*y, x*z, y*z
p(7:9) = 2*v(7:9);
p(10) = v(10);

% Get the explicit form of parameters
[center,radii,R] = ellipsoid_im2ex(p);

end

function [center,radii,R] = ellipsoid_im2ex(v)
% Cast ellipsoid defined with implicit parameter vector to explicit form.
% The implicit equation of a general ellipse is
% F(x,y,z) = Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz - 1 = 0
%
% Input arguments:
% v:
%    the 10 parameters describing the ellipsoid algebraically
% Output arguments:
% center:
%    ellispoid center coordinates [cx; cy; cz]
% ax:
%    ellipsoid semi-axes (radii) [a; b; c]
% R:
%    ellipsoid rotation (radii directions as rows of the 3x3 matrix)
%
% See also: ellipse_im2ex

% Copyright 2011 Levente Hunyadi

% eliminate times two from rotation and translation terms
v = v(:);
v(4:9) = 0.5*v(4:9);

% find the algebraic form of the ellipsoid (quadratic form matrix)
Q = [ v(1) v(4) v(5) v(7); ...
      v(4) v(2) v(6) v(8); ...
      v(5) v(6) v(3) v(9); ...
      v(7) v(8) v(9) v(10) ];

% find the center of the ellipsoid
center = Q(1:3,1:3) \ -v(7:9);

% form the corresponding translation matrix
T = eye(4,4);
T(4, 1:3) = center;

% translate to the center
S = T * Q * T';

% solve the eigenproblem
[evecs, evals] = eig( S(1:3,1:3) );
radii = real(sqrt( -S(4,4) ./ diag(evals) ));

% Rotation matrix
R = evecs';

end
