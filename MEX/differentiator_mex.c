#include <math.h>
#include "mex.h"

#define MIN(a,b) ((a) > (b) ? (b) : (a))
#define MAX(a,b) ((a) < (b) ? (b) : (a))
#define POW2(a) ((a)*(a))

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  double *x, *y, *dy, *coefs;
  double x_3, x_2, x_1, x0, x1, x2, x3, dx;
  double y_3, y_2, y_1, y0, y1, y2, y3;
  int npts, niter;
  bool is_super=false, *tmp_bool;
  double super_coefs[] = {2*58.0/252.0, 2*2*67.0/252.0, -2*3*22.0/252.0};
  double lanczos_coefs[] = {2*1.0/28.0, 2*2*2.0/28.0, 2*3*3.0/28.0};
  register int i, j;

  if (nrhs < 2) {
    mexErrMsgTxt("Not enough input arguments (2 is the minimum) !");
  }

  x = mxGetPr(prhs[0]);
  y = mxGetPr(prhs[1]);

  npts = mxGetM(prhs[0]);
  niter = mxGetN(prhs[0]);
  if (mxGetNumberOfElements(prhs[1]) != npts*niter) {
    mexErrMsgTxt("x and y should be matrices of the same size !");
  }

  if (nrhs > 2) {
    tmp_bool = mxGetLogicals(prhs[2]);
    is_super = tmp_bool[0];
  }

  if (is_super) {
    coefs = super_coefs;
  } else {
    coefs = lanczos_coefs;
  }

  plhs[0] = mxCreateDoubleMatrix(npts, niter, mxREAL);
  dy = mxGetPr(plhs[0]);

  for (i = 0; i < niter; i++) {
    x0 = x[0];
    x1 = x[1];
    x2 = x[2];
    x3 = x[3];

    dx = x1 - x0;

    x_3 = 2*x0 - x2 - dx;
    x_2 = 2*x0 - x1 - dx;
    x_1 = x0 - dx;

    y_3 = y[2];
    y_2 = y[1];
    y_1 = y[0];
    y0 = y[0];
    y1 = y[1];
    y2 = y[2];
    y3 = y[3];

    for (j = 0; j < npts-4; j++) {
      dy[j] = coefs[0] * (y1 - y_1) / (x1 - x_1) + coefs[1] * (y2 - y_2) / (x2 - x_2) + + coefs[2] * (y3 - y_3) / (x3 - x_3);

      x_3 = x_2;
      x_2 = x_1;
      x_1 = x0;
      x0 = x1;
      x1 = x2;
      x2 = x3;

      x3 = x[j+4];

      y_3 = y_2;
      y_2 = y_1;
      y_1 = y0;
      y0 = y1;
      y1 = y2;
      y2 = y3;

      y3 = y[j+4];
    }

    // npts-4
    dy[j] = coefs[0] * (y1 - y_1) / (x1 - x_1) + coefs[1] * (y2 - y_2) / (x2 - x_2) + + coefs[2] * (y3 - y_3) / (x3 - x_3);

    dx = x3 - x2;

    x_3 = x_2;
    x_2 = x_1;
    x_1 = x0;
    x0 = x1;
    x1 = x2;
    x2 = x3;

    x3 += dx;

    y_3 = y_2;
    y_2 = y_1;
    y_1 = y0;
    y0 = y1;
    y1 = y2;
    y2 = y3;

    //y3 = y[npts-1];

    j++;
    // npts-3
    dy[j] = coefs[0] * (y1 - y_1) / (x1 - x_1) + coefs[1] * (y2 - y_2) / (x2 - x_2) + + coefs[2] * (y3 - y_3) / (x3 - x_3);

    x_3 = x_2;
    x_2 = x_1;
    x_1 = x0;
    x0 = x1;
    x1 = x2;
    x2 = x3;

    x3 = 2*x[npts-1] - x[npts-2] + dx;

    y_3 = y_2;
    y_2 = y_1;
    y_1 = y0;
    y0 = y1;
    y1 = y2;
    y2 = y3;

    y3 = y[npts-2];

    j++;

    // npts-2
    dy[j] = coefs[0] * (y1 - y_1) / (x1 - x_1) + coefs[1] * (y2 - y_2) / (x2 - x_2) + + coefs[2] * (y3 - y_3) / (x3 - x_3);

    x_3 = x_2;
    x_2 = x_1;
    x_1 = x0;
    x0 = x1;
    x1 = x2;
    x2 = x3;

    x3 = 2*x[npts-1] - x[npts-3] + dx;

    y_3 = y_2;
    y_2 = y_1;
    y_1 = y0;
    y0 = y1;
    y1 = y2;
    y2 = y3;

    y3 = y[npts-3];

    j++;

    // npts-1
    dy[j] = coefs[0] * (y1 - y_1) / (x1 - x_1) + coefs[1] * (y2 - y_2) / (x2 - x_2) + + coefs[2] * (y3 - y_3) / (x3 - x_3);

    x += npts;
    y += npts;
    dy += npts;
  }

  return;
}
