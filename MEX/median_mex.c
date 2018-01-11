#include <math.h>
#include <string.h> 
#include "ctmf.h"
#include "mex.h"

#include "ctmf.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  int h, w, nelem, i, niter = 1, radius = 1; 
  unsigned char *tmp_img, *median_img, *tmp_ptr;
  double *img, value, mymin, mymax, scaling_factor;

  if (nrhs < 1) {
    mexErrMsgTxt("Not enough input arguments (1 is the minimum, 3 is the maximum) !");
  } else if (nrhs == 2) {
    radius = (int) mxGetScalar(prhs[1]);

    if (mxGetNumberOfElements(prhs[1]) > 1) {
      radius = radius / 2;
    }
  } else if (nrhs == 3) {
    radius = (int) mxGetScalar(prhs[1]);

    if (mxGetNumberOfElements(prhs[1]) > 1) {
      radius = radius / 2;
    }

    niter = (int) mxGetScalar(prhs[2]);
  }

  img = mxGetPr(prhs[0]);
  h = mxGetM(prhs[0]);
  w = mxGetN(prhs[0]);
  nelem = h*w;

  if ((median_img = mxCalloc(nelem, sizeof(unsigned char))) == NULL) {
    mexErrMsgTxt("Memory allocation failed !");
  }
  if ((tmp_img = mxCalloc(nelem, sizeof(unsigned char))) == NULL) {
    mexErrMsgTxt("Memory allocation failed !");
  }

  mymax = mymin = img[0];
  for (i = 1; i < nelem; i++) {
    if (img[i] < mymin) {
      mymin = img[i];
    } else if (img[i] > mymax) {
      mymax = img[i];
    }
  }
  scaling_factor = 255 / (mymax - mymin);

  for (i = 0; i < nelem; i++){
    if (mxIsNaN(img[i])) {
      median_img[i] = 0;
    } else {
      value = ceil(scaling_factor*(img[i] - mymin));

      median_img[i] = (unsigned char) value;
    }
  }

  for (i = 0; i < niter; i++) {
    tmp_ptr = tmp_img;
    tmp_img = median_img;
    median_img = tmp_ptr;

    ctmf(tmp_img, median_img, h, w, h, h, radius, 1, 3*1024*1024);
  }

  mxFree(tmp_img);
  
  plhs[0] = mxCreateDoubleMatrix(h, w, mxREAL);
  img = mxGetPr(plhs[0]);

  for (i=0;i < nelem; i++) {
    img[i] = ((double) ((double) median_img[i]) / scaling_factor) + mymin;
  }

  mxFree(median_img);

  return;
}
