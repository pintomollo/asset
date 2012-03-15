#include <math.h>
#include "mex.h"

static double correl(double *x, int pos, int l) {

  int i, n;
  double mx, my, sx, sy, sxy, xv, yv;

  n = ceil((double)l/2.0);

  /* Calculate the mean of the two series x[], y[] */
  mx = 0;
  my = 0;   
  for (i=0;i<n;i++) {
    mx += x[(l + pos - i) % l];
    my += x[(i + pos) % l];
  }
  mx /= n;
  my /= n;

  /* Calculate the denominator */
  sx = 0;
  sy = 0;
  sxy = 0;
  for (i=0;i<n;i++) {
    xv = x[(l + pos - i) % l] - mx;
    yv = x[(i + pos) % l] - my;
    sx += xv * xv;
    sy += yv * yv;
    sxy += xv * yv;
  }

  return -(sxy / sqrt(sx*sy));
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int h, w, i, j, *start, *length;
  double *img, *corr_img, nan_val = mxGetNaN();

  if (nrhs < 1) {
    mexErrMsgTxt("No input argument !");
  } else {

    img = mxGetPr(prhs[0]);
    h = mxGetM(prhs[0]);
    w = mxGetN(prhs[0]);

    plhs[0] = mxCreateDoubleMatrix(h, w, mxREAL);
    corr_img = mxGetPr(plhs[0]);

    if ((start = mxCalloc(w, sizeof(int))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((length = mxCalloc(w, sizeof(int))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }

    for (i = 0; i < w; i++) {
      for (j = 0; j < h; j++) {
        if (mxIsFinite(img[(i*h) + j])) {
          start[i] = j;

          break;
        } else {
          corr_img[(i*h) + j] = nan_val;
        }
      }

      for (j = h-1; j >= start[i]; j--) {
        if (mxIsFinite(img[(i*h) + j])) {
          length[i] = j - start[i] + 1;
          break;
        } else {
          corr_img[(i*h) + j] = nan_val;
        }
      }
      for (j = 0; j < length[i]; j++) {
        corr_img[(i*h) + j + start[i]] = correl(img + (i*h) + start[i], j, length[i]);
      }
    }

    mxFree(start);
    mxFree(length);
  }

  return;
}
