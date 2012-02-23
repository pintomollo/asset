#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  double *x_indx, *y_indx, *z_indx, *plane_indx, *errs, *center, *axes, orient, *res;
  double x, y, z, tmp, corient, sorient, a, b, z_coef, dist, accum;
  int page_size = 500, count, index=0, prev_plane, i;
  mwSize m, n;

  if (nrhs != 4) {
    mexErrMsgTxt("Ellipse distance needs 4 input arguments !");
  } else {
    m = mxGetM(prhs[0]); 
    n = mxGetN(prhs[0]);

    center = mxGetPr(prhs[1]);
    axes = mxGetPr(prhs[2]);
    orient = mxGetScalar(prhs[3]);

    if (n < 4) {
      mexErrMsgTxt("Ellipse coordinates should be a Nx4 matrix !");
    } else {
      x_indx = mxGetPr(prhs[0]);
      y_indx = x_indx + m;
      z_indx = y_indx + m;
      plane_indx = z_indx + m;
      prev_plane = plane_indx[0];
      corient = cos(orient);
      sorient = sin(orient);

      if ((errs = mxCalloc(page_size, sizeof(double))) == NULL) {
        mexErrMsgTxt("Memory allocation failed !");
      }

      for (i=0; i<m; i++){
        if (prev_plane != plane_indx[i]) {
          errs[index] = accum / (double)count;

          accum = 0;
          count = 0;
          index++;

          if ((index % page_size) == 0) {
            if ((errs = mxRealloc(errs, (index+page_size) * sizeof(double))) == NULL) {
              mexErrMsgTxt("Memory allocation failed !");
            }
          }
        }

        z = z_indx[i] - center[2];
        if (fabs(z) > axes[2]) {
          accum += 1;
        } else {
          z_coef = sqrt(1 - (pow(z,2) / pow(axes[2],2)));
          a = z_coef * axes[0];
          b = z_coef * axes[1];

          x = x_indx[i] - center[0];
          y = center[1] - y_indx[i];

          tmp = (x*corient + y*sorient);
          y = (y*corient - x*sorient);
          x = tmp;

          dist = fabs(hypot(x/a, y/b) - 1);
          accum += dist;
        }

        count++;
        prev_plane = plane_indx[i];
      }
      errs[index] = accum / (double)count;

      plhs[0] = mxCreateDoubleMatrix(index+1, 1, mxREAL);
      res = mxGetPr(plhs[0]);

      for (i=0; i<=index; i++){
        res[i] = errs[i];
      }

      mxFree(errs);
    }
  }
}
