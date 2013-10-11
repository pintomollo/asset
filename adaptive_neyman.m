function [pval, Tan, mhat] = adaptive_neyman(x1, x2, max_dim)
% ADAPTIVE_NEYMAN Computes the p-value of two groups of curves to be identical
%   using the adaptive Neyman statistics.
%
%   PVAL = ADAPTIVE_NEYMAN(Y1, Y2) computes PVAL between the sets of curves Y1 and Y2.
%   Both Y1 abd Y2 should be matrices were each column is a curve. Note that because
%   every curve must have the same number of data points, the end of the longest set
%   will be truncated. Note also that because FFT does not handle NaNs, the corresponding
%   data points will be removed from the analysis.
%
%   [PVAL, TAN, MHAT] = ADAPTIVE_NEYMAN(...) returns in addition TAN and MHAT, the two
%   statistical quantities defines in [1].
%
%   [...] = ADAPTIVE_NEYMAN(..., MAX_DIM) permits to set MAX_DIM to a defined value.
%   It basically defines the number of fourier components to take into consideration
%   for the statistical test (see [1]).
%
% Reference:
% [1] Fan, J. and Lin, S.K. (1998), Test of Significance when data are curves., 
%     Journal of American Statistical Association, 93, 1007-1021.
%
% Adapted from the S-code from Jianqing Fan:
% http://orfe.princeton.edu/~jqfan/publications-software.html

  L = min(size(x1, 1), size(x2, 1));
  x1 = x1(1:L, :);
  x2 = x2(1:L, :);

  goods = ~any(isnan([x1 x2]), 2);
  x1 = x1(goods, :);
  x2 = x2(goods, :);

  L = size(x1, 1);
  L = 2^nextpow2(L);

  X = fft(x1, L);
  Y = fft(x2, L);

  if (nargin == 2)
    max_dim = round(L/2);
    %max_dim = L;
  end
  max_dim = min(max_dim(1), L);

  if (max_dim <= 0)
    max_dim = 1;
  end

  R1 = real(X); I1 = imag(X);
  R2 = real(Y); I2 = imag(Y);

  X = R1; Y = R2;

  index = [2:2:L];
  X(index,:) = R1(2:floor(L/2+1),:);
  Y(index,:) = R2(2:floor(L/2+1),:);

  index = [3:2:L];
  X(index,:) = I1(2:floor(L/2+0.5),:);
  Y(index,:) = I2(2:floor(L/2+0.5),:);

  R1 = mean(X, 2);
  R2 = mean(Y, 2);
  I1 = var(X, 0, 2);
  I2 = var(Y, 0, 2);

  stddiff = (R1-R2)./sqrt(I1/size(X,2)+I2/size(Y,2));
  stddiff = stddiff(1:max_dim);
  stddiff = cumsum(stddiff.^2 - 1)./sqrt(2*[1:max_dim].');

  [stddiff, index] = sort(stddiff);
  Tan = sqrt(2* log(log(max_dim))) * stddiff - (2*log(log(max_dim)) + ...
        0.5*log(log(log(max_dim))) - 0.5*log(4*pi));

  % Original S-code
  %  mhat <- order(stddiff)[max_dim] 
  %  ANT <-  sqrt(2* log(log(max_dim))) *stddiff - (2*log(log(max_dim))+ 
  %      0.5*log(log(log(max_dim))) - 0.5*log(4*pi) )

  Tan = Tan(max_dim);
  mhat = index(max_dim);

  pval = adaptive_neyman_quantile_mex(max_dim, Tan);

  return;
end
