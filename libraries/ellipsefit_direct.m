function p = ellipsefit_direct(x,y)
% Direct least squares fitting of ellipses.
%
% Input arguments:
% x,y;
%    x and y coodinates of 2D points
%
% Output arguments:
% p:
%    a 6-parameter vector of the algebraic ellipse fit with
%    p(1)*x^2 + p(2)*x*y + p(3)*y^2 + p(4)*x + p(5)*y + p(6) = 0
%
% References:
% Andrew W. Fitzgibbon, Maurizio Pilu and Robert B. Fisher, "Direct Least
%    Squares Fitting of Ellipses", IEEE Trans. PAMI 21, 1999, pp476-480.

% Copyright 2011 Levente Hunyadi

narginchk(2,2);
validateattributes(x, {'numeric'}, {'real','nonempty','vector'});
validateattributes(y, {'numeric'}, {'real','nonempty','vector'});
x = x(:);
y = y(:);

% normalize data
mx = mean(x);
my = mean(y);
sx = (max(x)-min(x))/2;
sy = (max(y)-min(y))/2;
smax = max(sx,sy);
sx = smax;
sy = smax;
x = (x-mx)/sx;
y = (y-my)/sy;

% build design matrix
D = [ x.^2  x.*y  y.^2  x  y  ones(size(x)) ];

% build scatter matrix
S = D'*D;

% build 6x6 constraint matrix
C = zeros(6,6);
C(1,3) = -2;
C(2,2) = 1;
C(3,1) = -2;

% Do the actual fit
p = ellipsefit_robust(S,-C);

% unnormalize
p(:) = ...
[ p(1)*sy*sy ...
; p(2)*sx*sy ...
; p(3)*sx*sx ...
; -2*p(1)*sy*sy*mx - p(2)*sx*sy*my + p(4)*sx*sy*sy ...
; -p(2)*sx*sy*mx - 2*p(3)*sx*sx*my + p(5)*sx*sx*sy ...
; p(1)*sy*sy*mx*mx + p(2)*sx*sy*mx*my + p(3)*sx*sx*my*my - p(4)*sx*sy*sy*mx - p(5)*sx*sx*sy*my + p(6)*sx*sx*sy*sy ...
];

p = p ./ norm(p);

end

function p = ellipsefit_robust(R, Q)
% Constrained ellipse fit by solving a modified eigenvalue problem.
% The method is numerically stable.
%
% Input arguments:
% R:
%    positive semi-definite data covariance matrix
% Q:
%    constraint matrix in parameters x^2, xy, y^2, x, y and 1.
%
% Output arguments:
% p:
%    estimated parameters (taking constraints into account)

% References:
% Radim Halir and Jan Flusser, "Numerically stable direct least squares fitting of
%    ellipses", 1998

% Copyright 2012 Levente Hunyadi

S1 = R(1:3,1:3);     % quadratic part of the scatter matrix
S2 = R(1:3,4:6);     % combined part of the scatter matrix
S3 = R(4:6,4:6);     % linear part of the scatter matrix
T = -(S3 \ S2');     % for getting a2 from a1
M = S1 + S2 * T;     % reduced scatter matrix
M = Q(1:3,1:3) \ M;  % premultiply by inv(C1), e.g. M = [M(3,:)./2 ; -M(2,:) ; M(1,:)./2] for an ellipse
[evec,~] = eig(M);   % solve eigensystem

% evaluate a'*C*a, e.g. cond = 4 * evec(1,:).*evec(3,:) - evec(2,:).^2 for an ellipse
cond = zeros(1,size(evec,2));
for k = 1 : numel(cond)
    cond(k) = evec(:,k)'*Q(1:3,1:3)*evec(:,k);
end

% eigenvector for minimum positive eigenvalue
evec = evec(:,cond > 0);
cond = cond(cond > 0);
[~,ix] = min(cond);
p1 = evec(:,ix);  % eigenvector for minimum positive eigenvalue

% ellipse coefficients
p = [p1 ; T * p1];

end
