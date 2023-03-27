#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "gaussian_smooth.h"
#include "mex.h"

#include "gaussian_smooth.c"

/* A Matlab wrapper for the code of Anthony Gabrielson (see gaussian_smooth.c)*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  /* Declare a few variables. */
  mwSize h, w, c, dims[3];
  const mwSize *size;
  int j, offset_img;
  double *img, sigma;

  /* No flexibility here, we want both the image and sigma ! */
  if (nrhs < 2) {
    mexErrMsgTxt("Not enough input arguments (2 are required) !");
  } else if (!mxIsDouble(prhs[0])) {
    mexErrMsgTxt("Input array is not of type Double");
  }

  /* Get sigma. */
  sigma = mxGetScalar(prhs[1]);

  // The size of the image
  size = mxGetDimensions(prhs[0]);
  h = size[0];
  w = size[1];
  c = mxGetNumberOfElements(prhs[0]) / (h*w);

  // Prepare the output
  dims[0] = h;
  dims[1] = w;
  dims[2] = c;

  // Some temporary values
  offset_img = h*w;

  // Create the array
  plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);

  /* Create the output image on which we'll work directly. */
  img = mxGetPr(plhs[0]);

  /* Copy the input to the working image. */
  memcpy(img, mxGetPr(prhs[0]), h*w*c*sizeof(double)); 

  /* Verify that sigma is valid, and let's go ! */
  if (sigma <= 0) {
    mexWarnMsgTxt("Gaussian smoothing with invalid sigma !");
  } else {
    for (j=0; j < c; j++) {
      gaussian_smooth(img, w, h, sigma);
      img += offset_img;
    }
  }

  return;
}
