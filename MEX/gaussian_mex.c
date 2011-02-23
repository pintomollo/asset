#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include "gaussian_smooth.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  int h, w; 
  double *img, sigma;

  if (nrhs < 2) {
    mexErrMsgTxt("Not enough input arguments (2 are required) !");
  } else if (!mxIsDouble(prhs[0])) {
    mexErrMsgTxt("Input array is not of type Double");
  }
    
  sigma = mxGetScalar(prhs[1]);

  h = mxGetM(prhs[0]);
  w = mxGetN(prhs[0]);

  plhs[0] = mxCreateDoubleMatrix(h, w, mxREAL);
  img = mxGetPr(plhs[0]);
  memcpy(img, mxGetPr(prhs[0]), h*w*sizeof(double)); 

  if (sigma == 0) {
    mexWarnMsgTxt("Gaussian smoothing with simga equal to 0 !");
  } else {
    gaussian_smooth(img, w, h, sigma);
  }

  return;
}
