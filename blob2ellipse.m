function [center, axes_length, orientation] = blob2ellipse(pixels)
 %   Find the ellipse that has the same normalized second central moments as
 %   the region.  Compute the axes lengths, orientation, and eccentricity of
 %   the ellipse.  Ref: Haralick and Shapiro, Computer and Robot Vision vol
 %   I, Addison-Wesley 1992, Appendix A.
 
 if(size(pixels,2)~=2)
    inds = find(logical(pixels));
    [rows,cols] = ind2sub(size(pixels),inds);
    pixels = [cols rows];
 end
 
 center = mean(pixels,1)';

 % Assign X and Y variables so that we're measuring orientation
 % counterclockwise from the horizontal axis.
 
 xbar = center(1);
 ybar = center(2);
 
 x = pixels(:,1) - xbar;
 y = -(pixels(:,2) - ybar);
 % This is negative for the orientation calculation (measured
 % in the counter-clockwise direction).
 
 N = length(x);
 
 % Calculate normalized second central moments for the region.
 % 1/12 is the normalized second central moment of a pixel
 % with unit length.
 uxx = sum(x.^2)/N + 1/12;
 uyy = sum(y.^2)/N + 1/12;
 uxy = sum(x.*y)/N;
 
 % Calculate major axis length and minor axis length.
 common = sqrt((uxx - uyy)^2 + 4*uxy^2);
 axes_length = sqrt(2)*sqrt(uxx + uyy + common*[1; -1]);

 % Calculate orientation.
 if (uyy > uxx)
   num = uyy - uxx + sqrt((uyy - uxx)^2 + 4*uxy^2);
   den = 2*uxy;
 else
   num = 2*uxy;
   den = uxx - uyy + sqrt((uxx - uyy)^2 + 4*uxy^2);
 end
 if (num == 0) && (den == 0)
  orientation = 0;
 else
  orientation = atan2(num,den);
  orientation = 2*pi*(orientation<0) + orientation;
 end

 return;
end
