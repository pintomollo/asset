function Prob = chi2cdf(x,df)

% Purpose: 
% CDF of chi2 distribution
% -----------------------------------
% Density:
% f(x) = 1/(G(a) * b^a) * x^(a-1) * exp(-x/b), where a = df/2, b = 2
% E(X) = a*b, Var(X) = a*b^2
% -----------------------------------
% Algorithm: 
% Use incomplete gamma function 
% -----------------------------------
% Usage:
% x = points of evaluation
% df = degree of freedom parameter
% -----------------------------------
% Returns:
% Prob = Probability evaluated at points x
% -----------------------------------
% Notes:
%
% Written by Hang Qian, Iowa State University
% Contact me:  matlabist@gmail.com


Prob = gammainc(x./2,df/2);

