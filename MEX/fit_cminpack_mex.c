/* ////////////////////////////////////////////////////////////////////////////////
//
//  C code to optimize the front function used to quantify the cortical signal
//  in ASSET. This is a copy of the levmar.c function provided with the levmar
//  library, which I simplified for my very specific usage.
//
//  Simon Blanchoud
//  Naef & Gonczy labs
//  23.12.2011
//
//  Matlab MEX file for the Levenberg - Marquardt minimization algorithm
//  Copyright (C) 2007  Manolis Lourakis (lourakis at ics forth gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//////////////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>

#include <mex.h>
#include "cminpack_lmdif.h"

typedef struct {double *pos;
                double *values;
                double *smooth;
                double *lb;
                double *ub;
                double coef;
                } user_data;

static int func(void *d, int n, int m, double *p, double *hx, int iflag)
{
user_data *data = (user_data*)d;
int i, min_indx = 0, max_indx = n-1, params_inside=0;
double ampl, min_pos, max_pos, width, center, variance, base, slope;

  for (i=0;i<m;i++){
    params_inside += (int)(p[i] >= data->lb[i] && p[i] <= data->ub[i]);
  }

  if (params_inside==m) {
    ampl = p[0];
    center = p[1];
    width = data->coef*p[2];
    min_pos = center - width;
    max_pos = center + width;
    variance = -1/(2*p[2]*p[2]);

    memcpy(hx, data->smooth, n*sizeof(double));

    for (i = 0; i < n; i++) {
      hx[i] += ampl * exp((data->pos[i] - center)*(data->pos[i] - center)*variance) - data->values[i];
      if (data->pos[i] < min_pos) {
        min_indx++;
      }
      if (data->pos[i] > max_pos) {
        max_indx--;
      }
    }
    min_indx = min_indx > 1 ? min_indx-1 : min_indx;
    max_indx = max_indx < n-1 ? max_indx+1 : max_indx;

    base = data->smooth[min_indx];
    slope = (data->smooth[max_indx] - base) / (max_indx - min_indx);

    for (i = min_indx+1; i < max_indx; i++) {

      base += slope;
      hx[i] += base - data->smooth[i];
    }
  } else {
    memset(hx, ((m-params_inside)/m)*((unsigned char)DBL_MAX), n*sizeof(double));
  }

  return iflag;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
register int i, j;
int status=0, *iwa, lwa;
double *p, *p0, *ret, *x, *pos, *smoothed, stop_tol, *dwa, *fvec, coef, max_pos;
int m, n, niter, itmax, nbounds, tmp_n;
double *lb=NULL, *ub=NULL; 
user_data data;

  /* parse input args; start by checking their number */
  if((nrhs!=7))
    mexErrMsgTxt("cminpack: 7 input arguments required");

  /* the second required argument must be a real row or column vector */
  if(!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
    mexErrMsgTxt("cminpack: p0 must be a real vector.");

  p0=mxGetPr(prhs[0]);
  m = mxGetM(prhs[0]);
  niter = mxGetN(prhs[0]);

  /* copy input parameter vector to avoid destroying it, directly to the 2nd output argument */
  plhs[1]=mxCreateDoubleMatrix(m, niter, mxREAL);
  p=mxGetPr(plhs[1]);
  memcpy(p, p0, m*niter*sizeof(double));

  /** x **/
  /* the third required argument must be a real row or column vector */
  if(!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetN(prhs[1])!=niter)
    mexErrMsgTxt("cminpack: data must be a real matrix.");

  x=mxGetPr(prhs[1]);
  n=mxGetM(prhs[1]);

  /** coef **/
  /* the fourth required argument must be a scalar */
  if(!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)
    mexErrMsgTxt("cminpack: coef must be a scalar.");
  coef=mxGetScalar(prhs[2]);

  /** lb, ub **/
  /* check if the next two arguments are real row or column vectors */
  if(mxIsDouble(prhs[3]) && !mxIsComplex(prhs[3])) {
    if(mxIsDouble(prhs[4]) && !mxIsComplex(prhs[4])) {
      if(mxGetM(prhs[3])!=m)
        mexErrMsgTxt("cminpack: lb must have as many elements as p0");
      if(mxGetM(prhs[4])!=m)
        mexErrMsgTxt("cminpack: ub must have as many elements as p0");

      lb=mxGetPr(prhs[3]);
      ub=mxGetPr(prhs[4]);

      nbounds = mxGetN(prhs[3]);

      if (nbounds > 1 && nbounds != niter)
        mexErrMsgTxt("cminpack: ub and lb must have as many elements as p0");
    }
  }

  /** pos **/
  /* the 6th argument must be a real row or column vector */
  if(!mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) || mxGetNumberOfElements(prhs[5])!=n)
    mexErrMsgTxt("cminpack: pos must be a real vector.");
  pos=mxGetPr(prhs[5]);

  /** smoothed **/
  /* the last required argument must be a real row or column matrix */
  if(!mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) || mxGetM(prhs[6])!=n || mxGetN(prhs[6])!=niter)
    mexErrMsgTxt("cminpack: smoothed data must be a matrix similar to the data.");
  smoothed=mxGetPr(prhs[6]);

  data.pos = pos;
  data.smooth = smoothed;
  data.values = x;
  data.lb = lb;
  data.ub = ub;
  data.coef = coef;

  lwa = 2*(m * n + m * 5 + n);
  fvec = (double*)malloc(n*sizeof(double));
  iwa = (int*)malloc(m*sizeof(int));
  dwa = (double*)malloc(lwa*sizeof(double));

  stop_tol = sqrt(DBL_EPSILON);

  for (i = 0; i < niter; i++) {

    tmp_n = n-1;
    for (j = 0; j < n; j++) {
      if (mxIsNaN(data.smooth[j])) {
        tmp_n = j;
        break;
      }
    }
    max_pos = data.pos[tmp_n];

    if (tmp_n < m || (p[1] + coef*p[2]) >= max_pos) {
      //mexWarnMsgTxt("cminpack: not enough valid points to perform the optimization!");
      status = 1;
      p[2] = -1;
    } else {
      status=cminpack_lmdif(func, &data, tmp_n, m, p, fvec, stop_tol, iwa, dwa, lwa);

      for (j = 0; j < m; j++){
        p[j] = (p[j] < data.lb[j]) ? data.lb[j] : p[j];
        p[j] = (p[j] > data.ub[j]) ? data.ub[j] : p[j];
      }

      if ((p[1] + coef*p[2]) >= max_pos) {
        p[2] = -1;
      }
    }

    p+=m;
    data.values += n;
    data.smooth += n;

    if (nbounds > 1) {
      data.lb += m;
      data.ub += m;
    }

    if(status==0)
      mexWarnMsgTxt("cminpack: optimization returned with an error!");
  }

  free(iwa);
  free(dwa);
  free(fvec);

  /* copy back return flag */
  /** ret **/
  plhs[0]=mxCreateDoubleMatrix(1, 1, mxREAL);
  ret=mxGetPr(plhs[0]);
  ret[0]=(double)status;

  return;
}
