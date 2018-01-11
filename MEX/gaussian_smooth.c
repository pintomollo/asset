#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include "mex.h"
#include "gaussian_smooth.h"

/*******************************************************************************
* Adapted from MathWorks :
*
* File ID: #20899
* Canny Edge Detection
* by Anthony Gabrielson
*
* 29 Jul 2008 (Updated 30 Jul 2008)
*
* Code covered by the BSD License 
* http://www.mathworks.com/matlabcentral/fileexchange/20899-canny-edge-detection
* 
*******************************************************************************/
/*******************************************************************************
* PROCEDURE: gaussian_smooth
* PURPOSE: Blur an image with a gaussian filter.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void gaussian_smooth(double *image, int rows, int cols, double sigma)
{
   int r, c, rr, cc,     /* Counter variables. */
      windowsize,        /* Dimension of the gaussian kernel. */
      center;            /* Half of the windowsize. */
   double *tempim,        /* Buffer for separable filter gaussian smoothing. */
         *kernel,        /* A one dimensional gaussian kernel. */
         dot,            /* Dot product summing variable. */
         sum;            /* Sum of the kernel weights variable. */

   /****************************************************************************
   * Create a 1-dimensional gaussian smoothing kernel.
   ****************************************************************************/
   make_gaussian_kernel(sigma, &kernel, &windowsize);
   center = windowsize / 2;

   /****************************************************************************
   * Allocate a temporary buffer image and the smoothed image.
   ****************************************************************************/
   if((tempim = mxCalloc(rows*cols, sizeof(double))) == NULL){
      mexErrMsgTxt("Memory allocation failed for the buffer image !");
   }

   /****************************************************************************
   * Blur in the x - direction.
   ****************************************************************************/
   for(r=0;r<rows;r++){
      for(c=0;c<cols;c++){
         dot = 0.0;
         sum = 0.0;
         for(cc=(-center);cc<=center;cc++){
            if(((c+cc) >= 0) && ((c+cc) < cols)){
               dot += image[r*cols+(c+cc)] * kernel[center+cc];
               sum += kernel[center+cc];
            }
         }
         tempim[r*cols+c] = dot/sum;
      }
   }

   /****************************************************************************
   * Blur in the y - direction.
   ****************************************************************************/
   for(c=0;c<cols;c++){
      for(r=0;r<rows;r++){
         sum = 0.0;
         dot = 0.0;
         for(rr=(-center);rr<=center;rr++){
            if(((r+rr) >= 0) && ((r+rr) < rows)){
               dot += tempim[(r+rr)*cols+c] * kernel[center+rr];
               sum += kernel[center+rr];
            }
         }
         image[r*cols+c] = dot/sum;
      }
   }

   mxFree(tempim);
   mxFree(kernel);
}

/*******************************************************************************
* PROCEDURE: make_gaussian_kernel
* PURPOSE: Create a one dimensional gaussian kernel.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void make_gaussian_kernel(double sigma, double **kernel, int *windowsize)
{
   int i, center;
   double x, fx, sum=0.0, sqrttwopi;

    sqrttwopi = sqrt(M_PI * 2);

   *windowsize = 1 + 2 * ceil(2.5 * sigma);
   center = (*windowsize) / 2;

   if((*kernel = mxCalloc((*windowsize), sizeof(double))) == NULL){
      mexErrMsgTxt("Memory allocation failed for kernel !");
   }

   for(i=0;i<(*windowsize);i++){
      x = (double)(i - center);
      fx = exp(-0.5*x*x/(sigma*sigma)) / (sigma * sqrttwopi);
      (*kernel)[i] = fx;
      sum += fx;
   }

   for(i=0;i<(*windowsize);i++) (*kernel)[i] /= sum;
}

