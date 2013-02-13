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

//FILE *fid;

//static void func(double *p, double *hx, int m, int n)
static int func(void *d, int n, int m, double *p, double *hx, int iflag)
{
user_data *data = (user_data*)d;
bool params_inside = true;
int i, min_indx = 0, max_indx = n-1;
double ampl, min_pos, max_pos, width, center, variance, base, slope;

  for (i=0;i<m;i++){
    params_inside = (params_inside && p[i] >= data->lb[i] && p[i] <= data->ub[i]);
  //  fprintf(fid,"%e ", p[i]);
  }
  //fprintf(fid,"\n");

  if (params_inside) {
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
  //      hx[i] += data->smooth[i];
        min_indx++;
      }
      if (data->pos[i] > max_pos) {
  //      hx[i] += data->smooth[i];
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
    memset(hx, (unsigned char)DBL_MAX, n*sizeof(double));
  }

  return iflag;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
register int i, j;
int status=0, *iwa, lwa;
double *p, *p0, *ret, *x, *pos, *smoothed, stop_tol, *dwa, *fvec, coef;
int m, n, niter, itmax, nbounds;
//double opts[LM_OPTS_SZ]={LM_INIT_MU, LM_STOP_THRESH, LM_STOP_THRESH, LM_STOP_THRESH, LM_DIFF_DELTA};
double *lb=NULL, *ub=NULL; 
user_data data;

  /* parse input args; start by checking their number */
  if((nrhs!=8))
    mexErrMsgTxt("cminpack: 8 input arguments required");
    
  /* note that in order to accommodate optional args, prhs & nrhs are adjusted accordingly below */

  /* the second required argument must be a real row or column vector */
  if(!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
    mexErrMsgTxt("cminpack: p0 must be a real vector.");

  p0=mxGetPr(prhs[0]);
  m = mxGetM(prhs[0]);
  niter = mxGetN(prhs[0]);

  /* copy input parameter vector to avoid destroying it */
  plhs[1]=mxCreateDoubleMatrix(m, niter, mxREAL);
  p=mxGetPr(plhs[1]);
  memcpy(p, p0, m*niter*sizeof(double));

  /** x **/
  /* the third required argument must be a real row or column vector */
  //|| !(mxGetM(prhs[1])==1 || mxGetN(prhs[1])==1))
  if(!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetN(prhs[1])!=niter)
    mexErrMsgTxt("cminpack: data must be a real matrix.");

  x=mxGetPr(prhs[1]);
  n=mxGetM(prhs[1]);

  /** coef **/
  /* the fourth required argument must be a scalar */
  if(!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)
    mexErrMsgTxt("cminpack: coef must be a scalar.");
  coef=mxGetScalar(prhs[2]);

  /** stop_tol **/
  /* the fourth required argument must be a scalar */
  if(!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxGetM(prhs[3])!=1 || mxGetN(prhs[3])!=1)
    mexErrMsgTxt("cminpack: itmax must be a scalar.");
  stop_tol=mxGetScalar(prhs[3]);

  stop_tol = sqrt(DBL_EPSILON);
//  opts[1] = stop_tol;
//  opts[2] = stop_tol;
//  opts[3] = stop_tol;

  /* arguments below this point are optional and their presence depends
   * upon the minimization type determined above
   */
  /** lb, ub **/
    /* check if the next two arguments are real row or column vectors */
  if(mxIsDouble(prhs[4]) && !mxIsComplex(prhs[4])) {
  //&& (mxGetM(prhs[4])==1 || mxGetN(prhs[4])==1)){
    if(mxIsDouble(prhs[5]) && !mxIsComplex(prhs[5])) {
    //&& (mxGetM(prhs[5])==1 || mxGetN(prhs[5])==1)){
      //if((__MAX__(mxGetM(prhs[4]), mxGetN(prhs[4])))!=m)
      if(mxGetM(prhs[4])!=m)
        mexErrMsgTxt("cminpack: lb must have as many elements as p0");
      //if((__MAX__(mxGetM(prhs[5]), mxGetN(prhs[5])))!=m)
      if(mxGetM(prhs[5])!=m)
        mexErrMsgTxt("cminpack: ub must have as many elements as p0");

      lb=mxGetPr(prhs[4]);
      ub=mxGetPr(prhs[5]);

      nbounds = mxGetN(prhs[4]);

      if (nbounds > 1 && nbounds != niter)
        mexErrMsgTxt("cminpack: ub and lb must have as many elements as p0");
    }
  }

  /** pos **/
  /* the 7th argument must be a real row or column vector */
  if(!mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) || mxGetNumberOfElements(prhs[6])!=n)
    mexErrMsgTxt("cminpack: pos must be a real vector.");
  pos=mxGetPr(prhs[6]);

  /** smoothed **/
  /* the last required argument must be a real row or column matrix */
  if(!mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]) || mxGetM(prhs[7])!=n || mxGetN(prhs[7])!=niter)
    mexErrMsgTxt("cminpack: smoothed data must be a matrix similar to the data.");
  smoothed=mxGetPr(prhs[7]);

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

  //fid = fopen("sampling_test.txt", "a+");
  //fprintf(fid,"---------------CMINPACK-------------------\n");

  for (i = 0; i < niter; i++) {
  //for (i = 0; i < 1; i++) {

    //status=dlevmar_bc_dif(func, p, x, m, n, lb, ub, NULL, itmax, opts, NULL, NULL, NULL);
    //status=lmdif1((void*) func, n, m,(const double*) p, x, (double)stop_tol, iwa, dwa, lwa);
    status=cminpack_lmdif(func, &data, n, m, p, fvec, stop_tol, iwa, dwa, lwa);

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

  //fclose(fid);

  /* copy back return results */
  /** ret **/
  plhs[0]=mxCreateDoubleMatrix(1, 1, mxREAL);
  ret=mxGetPr(plhs[0]);
  ret[0]=(double)status;

  return;
}
