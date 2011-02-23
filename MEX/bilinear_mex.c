#include <math.h>
#include "mex.h"

/* Bilinear interpolation */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  int i, xf, yf, xc, yc;
  double dxf, dyf, dxc, dyc, x, y;
  mwSize w, h, m, n, nvals;
  double *x_indx, *y_indx, *tmp, *img, *values;
  bool free_memory = false, circular_x = false, circular_y = false, *tmp_bool;

  if (nrhs < 2) {
    mexErrMsgTxt("Not enough input arguments (2 is the minimum) !");
  } else if (nrhs == 2) {
    m = mxGetM(prhs[1]); 
    n = mxGetN(prhs[1]);

    if (n == 2) {
      x_indx = mxGetPr(prhs[1]);
      y_indx = x_indx + m;
      
      n = 1;
    } else if (m == 2) {
      if ((x_indx = mxCalloc(n, sizeof(double))) == NULL) {
        mexErrMsgTxt("Memory allocation failed !");
      }
      if ((y_indx = mxCalloc(n, sizeof(double))) == NULL) {
        mexErrMsgTxt("Memory allocation failed !");
      }

      tmp = mxGetPr(prhs[1]);

      for (i = 0; i < n; i++) {
        x_indx[i] = tmp[i*2];
        y_indx[i] = tmp[i*2 + 1];
      }

      m = 1;

      free_memory = true;

    } else {
      mexErrMsgTxt("Not enough input arguments if both X and Y indexes are not included !");
    }
  } else if (nrhs == 3) {
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != mxGetM(prhs[2])*mxGetN(prhs[2])){
      mexErrMsgTxt("Both index vectors must have the same length");
    }
    x_indx = mxGetPr(prhs[1]);
    y_indx = mxGetPr(prhs[2]);

    m = mxGetM(prhs[1]); 
    n = mxGetN(prhs[1]);

  } else if (nrhs == 4) {
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != mxGetM(prhs[2])*mxGetN(prhs[2])){
      mexErrMsgTxt("Both index vectors must have the same length");
    }
    x_indx = mxGetPr(prhs[1]);
    y_indx = mxGetPr(prhs[2]);

    m = mxGetM(prhs[1]); 
    n = mxGetN(prhs[1]);

    tmp_bool = mxGetLogicals(prhs[3]);
    circular_x = tmp_bool[0];

    if (mxGetNumberOfElements(prhs[3]) > 1) {
      circular_y = tmp_bool[1];
    } else {
      circular_y = tmp_bool[0];
    }
  }
   
  plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  values = mxGetPr(plhs[0]);

  h = mxGetM(prhs[0]);
  w = mxGetN(prhs[0]);
  img = mxGetPr(prhs[0]);

  /* Bilinear interpolation*/
  for (i=0; i < m*n; i++) {
    x = x_indx[i] - 1;
    y = y_indx[i] - 1;

    xf = floor(x);
    yf = floor(y);

    dxf = x - xf;
    dyf = y - yf;

    if (dxf == 0) {
      xc = xf;
      dxc = 1;
    } else {
      xc = xf + 1;
      dxc = xc - x; 
    }

    if (dyf == 0) {
      yc = yf;
      dyc = 1;
    } else {
      yc = yf + 1;
      dyc = yc - y;
    }

    if (circular_x) {
      xf = fmod(xf, w);
      xc = fmod(xc, w);

      if (xf == 0) {
        xf = w-1;
      }
      if (xc == 0) {
        xc = w-1;
      }
    }

    if (circular_y) {
      yf = fmod(yf,h);
      yc = fmod(yc,h);

      if (yf == 0) {
        yf = h-1;
      }
      if (yc == 0) {
        yc = h-1;
      }
    }

    if (xf >= w || yf >= h || xc < 0 || yc < 0 || xc >= w || xf < 0 || yc >= h || yf < 0) {
      values[i] = mxGetNaN();

    } else {
      values[i] = img[xf*h + yf] * dxc * dyc + 
                  img[xc*h + yf] * dxf * dyc +
                  img[xf*h + yc] * dxc * dyf +
                  img[xc*h + yc] * dxf * dyf;
    } 
  }

  if (free_memory) {
    mxFree(x_indx);
    mxFree(y_indx);
  }

  return;
}
