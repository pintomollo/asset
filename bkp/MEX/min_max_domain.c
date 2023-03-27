#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  int i, j, p, pos, nnoise = 0;
  mwSize w, h, planes, half, ndims;
  const mwSize *sizes;
  double *img, *path, *noises, noise = 0, *column, ref_val, tmp_val;

  if (nrhs < 2) {
    mexErrMsgTxt("Not enough input arguments (2 is the minimum) !");
  } else if (nrhs == 3) {
    noises = mxGetPr(prhs[2]);
    nnoise = mxGetNumberOfElements(prhs[2]);
  }
  
  ndims = mxGetNumberOfDimensions(prhs[0]);
  sizes = mxGetDimensions(prhs[0]);
  plhs[0] = mxCreateNumericArray(ndims, sizes, mxDOUBLE_CLASS, mxREAL);

  h = sizes[0];
  w = sizes[1];
  planes = ndims > 2 ? sizes[2] : 1;

  img = mxGetPr(plhs[0]);
  memcpy(img, mxGetPr(prhs[0]), h*w*planes*sizeof(double)); 

  half = h/2;

  path = mxGetPr(prhs[1]);

  column = img;

  for (p=0; p<planes; p++) {
    if (nnoise == planes) {
      noise = noises[p];
    }
    for (i=0; i < w; i++) {
      if (mxIsNaN(path[i])) {
        pos = half - 1;
      } else {
        pos = half - (int)path[i] - 1;
      }
      pos = pos < 0 ? 0 : pos;
      pos = pos >= half-1 ? half-1 : pos;

      ref_val = column[pos];
      for (j=pos-1;j>=0;j--) {
        tmp_val = column[j];
        if (mxIsFinite(tmp_val)) {
          ref_val = ref_val > tmp_val + noise ? tmp_val : ref_val;
          column[j] = tmp_val > ref_val + noise ? ref_val : tmp_val;
        } else {
          column[j] = ref_val;
        }
      }
      /*
      ref_val = column[pos];
      for (j=pos+1;j<half;j++) {
        tmp_val = column[j];
        if (mxIsFinite(tmp_val)) {
          ref_val = ref_val < tmp_val - noise ? tmp_val : ref_val;
          column[j] = tmp_val < ref_val - noise ? ref_val : tmp_val;
        } else {
          column[j] = ref_val;
        }
      }*/

      column+=half;
      ref_val = column[pos];
      for (j=pos-1;j>=0;j--) {
        tmp_val = column[j];
        if (mxIsFinite(tmp_val)) {
          ref_val = ref_val > tmp_val + noise ? tmp_val : ref_val;
          column[j] = tmp_val > ref_val + noise ? ref_val : tmp_val;
        } else {
          column[j] = ref_val;
        }
      }
      /*
      ref_val = column[pos];
      for (j=pos+1;j<half;j++) {
        tmp_val = column[j];
        if (mxIsFinite(tmp_val)) {
          ref_val = ref_val < tmp_val - noise ? tmp_val : ref_val;
          column[j] = tmp_val < ref_val - noise ? ref_val : tmp_val;
        } else {
          column[j] = ref_val;
        }
      }*/
      column+=half;
    }
  }

  return;
}
