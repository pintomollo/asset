#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include "gaussian_smooth.h"
#include "mex.h"

#ifndef MAX
#define MAX(a,b) ((a) < (b) ? (b) : (a))
#endif

#ifndef EDGE_DIR
#define EDGE_DIR(a,b,c,d) \
    ((a) < (b) ? ((a) < (c) ? ((a) < (d) ? 1 : 4) : ((c) < (d) ? 3 : 4)) : \
                 ((b) < (c) ? ((b) < (d) ? 2 : 4) : ((c) < (d) ? 3 : 4)));
#endif


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  int h, w, i, j;
  double thresh = 0, *img, *tmp_img, *direct;
  double hedge, vedge, pdedge, ndedge, min_val = 5, max_val = 0;
  double angle, opp_angle, pi2 = M_PI/2, pi4 = M_PI/4, twopi = 2*M_PI;
  unsigned short int *dir;
  bool single = true;

  if (nrhs == 2) {
    thresh = mxGetScalar(prhs[1]);
  } else if (nrhs == 3) {
    thresh = mxGetScalar(prhs[1]);
    single = (bool) mxGetScalar(prhs[2]);
  }

  h = mxGetM(prhs[0]);
  w = mxGetN(prhs[0]);

  plhs[0] = mxCreateDoubleMatrix(h, w, mxREAL);
  img = mxGetPr(plhs[0]);
  memcpy(img, mxGetPr(prhs[0]), h*w*sizeof(double)); 

  gaussian_smooth(img, w, h, 1.5);

  if ((tmp_img = mxCalloc(h*w, sizeof(double))) == NULL) {
    mexErrMsgTxt("Memory allocation failed !");
  }
  if ((dir = mxCalloc(h*w, sizeof(unsigned short int))) == NULL) {
    mexErrMsgTxt("Memory allocation failed !");
  }
  if (nlhs > 1) {
    plhs[1] = mxCreateDoubleMatrix(h, w, mxREAL);
    direct = mxGetPr(plhs[1]);
  }  

  for (i=0; i < w; i++) {
    for (j=0; j < h; j++) {

      if (i > 2 && j < h-2 && i < w-2 && j> 2) {
        ndedge = fabs(img[(j+2) + h*(i-2)] + img[(j+1) + h*(i-1)]
                    - img[(j-1) + h*(i+1)] - img[(j-2) + h*(i+2)]);
      } else {
        ndedge = 0;
      }

      if (j > 2 && j < h-2) {
        vedge = fabs(img[(j-2) + h*i] + img[(j-1) + h*i]
                   - img[(j+1) + h*i] - img[(j+2) + h*i]);
      } else {
        vedge = 0;
      }

      if (i > 2 && j > 2 && i < w-2 && j < h-2) {
        pdedge = fabs(img[(j-2) + h*(i-2)] + img[(j-1) + h*(i-1)]
                    - img[(j+1) + h*(i+1)] - img[(j+2) + h*(i+2)]);
      } else {
        pdedge = 0;
      }

      if (i > 2 && i <w-2) {
        hedge = fabs(-img[j + h*(i-2)] - img[j + h*(i-1)]
                    + img[j + h*(i+1)] + img[j + h*(i+2)]);
      } else {
        hedge = 0;
      }

      tmp_img[j + h*i] = MAX(MAX(ndedge,vedge), MAX(pdedge, hedge)) / 2;
      dir[j + h*i] = EDGE_DIR(ndedge, vedge, pdedge, hedge);

      if (tmp_img[j + h*i] < min_val) min_val = tmp_img[j + h*i];
      if (tmp_img[j + h*i] > max_val) max_val = tmp_img[j + h*i];

      if (nlhs > 1) {
        angle = atan2(hedge, vedge);
        opp_angle = atan2(ndedge, pdedge) - pi4;
        if (angle > pi2 && opp_angle < -pi2) opp_angle += twopi; 

        direct[j + h*i] = (angle + opp_angle) / 2;
      }
    }
  }

  memcpy(img, tmp_img, h*w*sizeof(double));

  if (single) {
    for (i=0; i < w; i++) {
      for (j=0; j < h; j++) {
        switch (dir[j + h*i]) {
          case 1:
            if (i > 0 && j > 0 && tmp_img[j + h*i] < tmp_img[(j-1) + h*(i-1)]) {
              img[j + h*i] = 0;
              break;
            }
            if (i < w-1 && j < h-1 && tmp_img[j + h*i] <= tmp_img[(j+1) + h*(i+1)]) {
              img[j + h*i] = 0;
            }

            break;
          case 2:
            if (i > 0 && tmp_img[j + h*i] < tmp_img[j + h*(i-1)]) {
              img[j + h*i] = 0;
              break;
            }
            if (i < w-1 && tmp_img[j + h*i] <= tmp_img[j + h*(i+1)]) {
              img[j + h*i] = 0;
            }

            break;
          case 3:
            if (i > 0 && j < h-1 && tmp_img[j + h*i] < tmp_img[(j+1) + h*(i-1)]) {
              img[j + h*i] = 0;
              break;
            }
            if (i < w-1 && j > 0 && tmp_img[j + h*i] <= tmp_img[(j-1) + h*(i+1)]) {
              img[j + h*i] = 0;
            }

            break;
          case 4:
            if (j > 0 && tmp_img[j + h*i] < tmp_img[(j-1) + h*i]) {
              img[j + h*i] = 0;
              break;
            }
            if (j < h-1 && tmp_img[j + h*i] <= tmp_img[(j+1) + h*i]) {
              img[j + h*i] = 0;
            }

            break;
        }

        if (img[j + h*i] > 0) {
          img[j + h*i] = (img[j + h*i] - min_val) / (max_val - min_val);
          if (img[j + h*i] < thresh) img[j + h*i] = 0;
        } 
        if (nlhs > 1 && img[j + h*i] == 0) {
          direct[j + h*i] = 0;
        }
      }
    }
  } else {
    for (i = 0; i < h*w; i++) {
      img[i] = (img[i] - min_val) / (max_val - min_val);
      if (img[i] < thresh) {
        img[i] = 0;

        if (nlhs > 1) direct[i] = 0;
      }
    }
  }

  mxFree(tmp_img);
  mxFree(dir);

  return;
}
