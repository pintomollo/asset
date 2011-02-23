#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  double *values, *candidates, *datas, *prev_dist, *prev_dir; 
  double *bests, *emission, *trans;
  double alpha, beta, gamma, inv_alpha, inv_beta, inv_gamma, smooth, norm, twonorm;
  double real_dir, tmp_dir, res, min_dist, myInf, tmp_val;
  int nhood, npts, half, i, p, *indxs, index, min_indx;
  bool no_values;

  if (nrhs < 6) {
    mexErrMsgTxt("Scoring function requires 6 inputs");
  }

  values = mxGetPr(prhs[0]);
  candidates = mxGetPr(prhs[1]);
  datas = mxGetPr(prhs[2]);
  prev_dist = mxGetPr(prhs[3]);
  prev_dir = mxGetPr(prhs[4]);

  nhood = (int) mxGetScalar(mxGetField(prhs[5], 0, "nhood"));
  alpha = mxGetScalar(mxGetField(prhs[5], 0, "alpha"));
  beta = mxGetScalar(mxGetField(prhs[5], 0, "beta"));
  gamma = mxGetScalar(mxGetField(prhs[5], 0, "gamma"));
  inv_alpha = 1 - alpha;
  inv_beta = 1 - beta;
  inv_gamma = 1 - gamma;
  myInf = mxGetInf();

  no_values = (mxGetNumberOfElements(prhs[0]) == 0);
  npts = mxGetNumberOfElements(prhs[1]);
  half = nhood / 2;
  if (half != 0) {
    norm = 1/((double)half);
    twonorm = norm / 2;
  } else {
    twonorm = 1;
    norm = 1;
  }

  plhs[0] = mxCreateDoubleMatrix(1, npts, mxREAL);
  bests = mxGetPr(plhs[0]);
  plhs[1] = mxCreateNumericMatrix(1, npts, mxUINT32_CLASS, mxREAL);
  indxs = (int *) mxGetData(plhs[1]);

  if (nlhs > 2) {
    plhs[2] = mxCreateDoubleMatrix(nhood, npts, mxREAL);
    emission = mxGetPr(plhs[2]);
  }
  if (nlhs > 3) {
    if (no_values) {
      plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
      trans = mxGetPr(plhs[3]);
      trans[0] = 1;
    } else {
      plhs[3] = mxCreateDoubleMatrix(1, nhood, mxREAL);
      trans = mxGetPr(plhs[3]);
    }
  }

  for (p = 0; p < npts; p++) {
    if (mxIsNaN(candidates[p]) || mxIsNaN(datas[p])) {
      bests[p] = myInf;
      indxs[p] = p + 1;

    } else {

      min_dist = myInf;
      min_indx = p + 1;

      for (i = -half; i <= half; i++) {
        if (nlhs > 3 && p == 0 && !no_values) {
          trans[i+half] = exp(-fabs(i)*norm*alpha*beta*gamma);
        }

        index = fmod(npts+p+i,npts);

        if ((no_values || !mxIsNaN(values[index])) && prev_dist[index] != myInf) {

          if (no_values) {
            smooth = 0;
          } else {
            real_dir = prev_dir[index] - index - 1;
            tmp_dir = prev_dir[index] - index - npts - 1;
            if (fabs(real_dir) > fabs(tmp_dir)) {
              real_dir = tmp_dir;
            }

            smooth = (fabs(i)*norm*gamma + 
                      (fabs(i - real_dir) * twonorm) * inv_gamma)
                       *beta
                    + fabs(candidates[p] - values[index]) * inv_beta;
          }

          res = alpha*smooth + inv_alpha*datas[p] + prev_dist[index];
          if (res < min_dist) {
            min_dist = res;
            min_indx = index + 1;
          }

          if (nlhs > 2) {
            if (no_values) {
              emission[i+half + p*nhood] = exp(-inv_alpha * datas[p]);
            } else {
              emission[i+half + p*nhood] = exp(-alpha * 
                        ((fabs(i - real_dir) * twonorm) 
                        * inv_gamma * gamma + 
                        fabs(candidates[p] - values[index]) * inv_beta) 
                        - inv_alpha * datas[p]);
            }
          }
        }
      }

      bests[p] = min_dist;
      indxs[p] = min_indx;
    }
  }

  return;
}
