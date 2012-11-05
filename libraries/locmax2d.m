function fImg=locmax2d(img,mask,keepFlat)
%LOCALMAX searches for local maxima in an image
%
%    SYNOPSIS fImg=(img,mask)
%
%    INPUT    img    image matrix
%             mask   [m n] defines the operator window dimensions
%             keepFlat Optional input variable to choose whether to remove
%                      "flat" maxima or to keep them. Default is 0, to
%                      remove them. - KJ
%
%    OUTPUT   fImg   map with all the local maxima (zeros elsewhere);
%                    the non-zero values contain the original value 
%                    of the image at that place
%
% NOTE: convert fImg to uint8 or uint16 for optimal display!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PARAMETER CHECK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
   error('Please define all parameters');
end

if nargin < 3 || isempty(keepFlat)
    keepFlat = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DEFINITIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make sure the mask elements are odd numbers (only then, the local max operator is 
% properly defined)
indx = find(~mod(mask,2));
mask(indx) = mask(indx) + 1;

% apply a max filter
fImg = ordfilt2(img,prod(mask),ones(mask));
if keepFlat == 0 %change made by KJ
    fImg2 = ordfilt2(img,prod(mask)-1,ones(mask));
    fImg(fImg2==fImg)=0;
end

% take only those positions where the max filter and the original image value
% are equal -> this is a local maximum
fImg(fImg ~= img) = 0;

% set image border to zero
auxM = zeros(size(img));
auxM(fix(mask(1)/2)+1:end-fix(mask(1)/2),fix(mask(2)/2)+1:end-fix(mask(2)/2)) = ...
    fImg(fix(mask(1)/2)+1:end-fix(mask(1)/2),fix(mask(2)/2)+1:end-fix(mask(2)/2));
fImg=auxM;
