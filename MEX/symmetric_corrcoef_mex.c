#include <math.h>
#include "mex.h"

#ifndef MIN
#define MIN(a,b) ((a) > (b) ? (b) : (a))
#endif

static double correl(double *x, int pos, int n, int l) {

  int i, n_max;
  double mx, my, sx, sy, sxy, xv, yv, corr;

  n_max = ceil((double)l/2.0);
  n = MIN(n, n_max);

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

  if (sxy == 0) {
    corr = 0;
  } else {
    corr = -(sxy / sqrt(sx*sy));
  }

  return corr;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int h, w, i, j, *start, *length;
  double *img, *lengths, *corr_img, nan_val = mxGetNaN(), window_size;
  bool *edges, multi_scales = false;

  if (nrhs < 1) {
    mexErrMsgTxt("No input argument !");
  } else {

    img = mxGetPr(prhs[0]);
    h = mxGetM(prhs[0]);
    w = mxGetN(prhs[0]);

    if (nrhs == 2) {
      if (mxIsLogical(prhs[1])) {
        mexErrMsgTxt("Not implemented yet !");
      } else {
        if (mxGetNumberOfElements(prhs[1]) == 1) {
          window_size = mxGetScalar(prhs[1]);
        } else if (mxGetNumberOfElements(prhs[1]) == 0) {
          window_size = ceil((double)w/2.0);
        } else {
          lengths = mxGetPr(prhs[1]);
          multi_scales = true;
        }
      }
    } else {
      window_size = ceil((double)w/2.0);
    }

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
        //corr_img[(i*h) + j + start[i]] = correl(img + (i*h) + start[i], j, length[i]);
        if (multi_scales) {
          corr_img[(i*h) + j + start[i]] = correl(img + (i*h) + start[i], j, (int)lengths[(i*h) + j + start[i]], length[i]);
        } else {
          corr_img[(i*h) + j + start[i]] = correl(img + (i*h) + start[i], j, window_size, length[i]);
        }
      }
    }

    mxFree(start);
    mxFree(length);
  }

  return;
}
