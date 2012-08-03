#include <math.h>
#include "mex.h"

/*
 * percentile.c - Selecting an item at a percentile in a list
 *
 * Copyright (C) 2008 by Boulder Laboratory for 3-Dimensional Electron
 * Microscopy of Cells ("BL3DEMC") and the Regents of the University of 
 * Colorado.  See dist/COPYRIGHT for full notice.
 *
 * Selects item number [s] (numbered from 1) out of [num] items in the array
 * [r], where items are considered in order from low to high.  [r] is partially
 * rearranged while finding the item.
 */
static double percentile(int s, double *r, int num)
{
  int lo = 0;
  int up = num - 1;
  int i, j;
  double temp;
  s--;
  while (up >= s && s >= lo) {
    i = lo;
    j = up;
    temp = r[s];
    r[s] = r[lo];
    r[lo] = temp;
    while (i < j) {
      while (r[j] > temp)
        j--;
      r[i] = r[j];
      while (i < j && r[i] <= temp)
        i++;
      r[j] = r[i];
    }
    r[i] = temp;
    if (s < i)
      up = i - 1;
    else
      lo = i + 1;
  }

  return r[s];
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  double *values, *candidates, *datas, *prev_dist, *prev_dir; 
  double *bests, *emission, *trans;
  double alpha, beta, gamma, inv_alpha, inv_beta, inv_gamma, smooth, norm, twonorm;
  double real_dir, tmp_dir, res, min_dist, myInf, tmp_val, spawn = 1;
  int nhood, npts, half, i, p, *indxs, index, min_indx, nvalids, spawn_index;
  bool no_values, transitions_done = false, is_circular = false, does_spawn = false;

  if (nrhs < 6) {
    mexErrMsgTxt("Scoring function requires at least 6 inputs");
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

  if (nrhs == 7) {
    if (mxIsLogical(prhs[6])) {
      is_circular = mxGetLogicals(prhs[6])[0];
    } else {
      spawn = mxGetScalar(prhs[6]);
    }
  } else if (nrhs == 8) {
    if (mxIsLogical(prhs[6])) {
      is_circular = mxGetLogicals(prhs[6])[0];
      spawn = mxGetScalar(prhs[7]);
    } else {
      spawn = mxGetScalar(prhs[6]);
      is_circular = mxGetLogicals(prhs[7])[0];
    }
  }

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

  if (spawn < 1 && spawn > 0) {
    does_spawn = true;
    nvalids = npts;
    for (p = 0; p < npts; p++) {
      if (prev_dist[p] == myInf) {
        nvalids--;
      }
    }

    if (myInf > 0) {
      spawn_index = (int)((spawn * (float)nvalids) + 1);
    } else {
      spawn_index = (int)((spawn * (float)nvalids) + 1 + npts - nvalids);
    }

    spawn = percentile(spawn_index, prev_dist, npts);
  }

  for (p = 0; p < npts; p++) {
    if (mxIsNaN(candidates[p]) || mxIsNaN(datas[p])) {
      bests[p] = myInf;
      indxs[p] = p + 1;

    } else {

      if (does_spawn) {
        min_dist = inv_alpha*datas[p] + spawn;
        min_indx = 0;
      } else {
        min_dist = myInf;
        min_indx = p + 1;
      }

      for (i = -half; i <= half; i++) {
        if (nlhs > 3 && !transitions_done && !no_values) {
          trans[i+half] = exp(-fabs(i)*norm*alpha*beta*gamma);
        }

        if (is_circular) {
          index = fmod(npts+p+i,npts);
        } else {
          index = p+i;
        }

        if (index >= 0 && index < npts) {
          if ((no_values || !mxIsNaN(values[index])) && prev_dist[index] != myInf) {
            if (no_values || (does_spawn && prev_dir[index] == 0)) {
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
      }
      if (nlhs > 3 && !transitions_done && !no_values) {
        transitions_done = true;
      }

      bests[p] = min_dist;
      indxs[p] = min_indx;
    }
  }

  return;
}
