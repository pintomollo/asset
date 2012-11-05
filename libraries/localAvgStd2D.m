% Computes the local average and standard deviation within a square window of side 'w'.

% Francois Aguet
% Last modified on 10/14/2010

%% Taken from spotDetector
% [frameInfo imgDenoised] = detectSpotsWT(img, S, dthreshold, postProcLevel)
%
% Performs detection of local intensity clusters through a combination of 
% multiscale products and denoising by iterative filtering from
% significant coefficients:
% Olivo-Marin, "Extraction of spots in biological images using multiscale products," Pattern Recoginition 35, pp. 1989-1996, 2002.
% Starck et al., "Image Processing and Data Analysis," Section 2.3.4, p. 73
%
% INPUTS:   img             : input image (2D array)
%           {S}             : postprocessing level.
%           {dthreshold}    : minimum allowed distance of secondary maxima in large clusters
%           {postProcLevel} : morphological post processing level for mask 

% Parts of this function are based on code by Henry Jaqaman.
% Francois Aguet, March 2010

function [avg sigma] = localAvgStd2D(image, w)

if mod(w+1, 2)
    error('The window length w should be an odd integer.');
end;

b = (w-1)/2;
image = padarray(image, [b b], 'replicate');

h = ones(1,w);
E = conv2(h/w, h, image, 'valid');
E2 = conv2(h, h, image.^2, 'valid');

sigma = E2 - E.^2;
sigma(sigma<0) = 0;
sigma = sqrt(sigma/(w*w - 1));
avg = E/w;
