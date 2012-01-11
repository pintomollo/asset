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

#include <levmar.h>

#include <mex.h>

#define __MAX__(A, B)     ((A)>=(B)? (A) : (B))
#define MIN_CONSTRAINED_BC    1

static void func(double *p, double *hx, int m, int n, void *adata)
{
int i;
double *pos=(double *)adata;
double min_pos, max_pos, width, center, variance, base, knot1, knot2, slope, mid_point;

  center = p[0];
  width = 1.79*p[1];
  min_pos = center - width;
  max_pos = center + width;
  variance = 2*pow(p[1], 2);

  /* Invagination */
  if (m == 6) {

    knot1 = p[3]*width + p[4];
    knot2 = -p[2]*width + p[4];
    mid_point = (knot1 + knot2) / 2;
    slope = (knot2 - knot1) / (2*width);

    for(i=0; i<n; ++i) {
      if (pos[i] < min_pos) {
        base = -p[3]*(pos[i] - center) + p[4];
      } else if (pos[i] < max_pos) {
        base = slope * (pos[i] - center) + mid_point;
      } else {
        base = -p[2]*(pos[i] - center) + p[4];
      }

      hx[i] = base + p[5] * exp(-pow(pos[i] - center, 2) / variance);
    }

  /* Front */
  } else {

    knot1 = p[4]*width + p[5];
    knot2 = p[5] - p[3]*(1 - exp(-p[2]*width));
    mid_point = (knot1 + knot2) / 2;
    slope = (knot2 - knot1) / (2*width);

    for(i=0; i<n; ++i) {
      if (pos[i] < min_pos) {
        base = -p[4]*(pos[i] - center) + p[5];
      } else if (pos[i] < max_pos) {
        base = slope * (pos[i] - center) + mid_point;
      } else {
        base = p[5] - p[3]*(1 - exp(-p[2]*(pos[i] - center)));
      }

      hx[i] = base + p[6] * exp(-pow(pos[i] - center, 2) / variance);
    }
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *Prhs[])
{
register int i;
register double *pdbl;
mxArray **prhs=(mxArray **)&Prhs[0], **ptrs;
int status, isrow_p0;
double *p, *p0, *ret, *x, *pos, stop_tol;
int m, n, itmax, mintype;
double opts[LM_OPTS_SZ]={LM_INIT_MU, LM_STOP_THRESH, LM_STOP_THRESH, LM_STOP_THRESH, LM_DIFF_DELTA};
double *lb=NULL, *ub=NULL; 

  /* parse input args; start by checking their number */
  if((nrhs!=7))
    mexErrMsgTxt("fit_front: 7 input arguments required");
    
  /* note that in order to accommodate optional args, prhs & nrhs are adjusted accordingly below */

  /* the second required argument must be a real row or column vector */
  if(!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || !(mxGetM(prhs[0])==1 || mxGetN(prhs[0])==1))
    mexErrMsgTxt("fit_front: p0 must be a real vector.");
  p0=mxGetPr(prhs[0]);

  /* determine if we have a row or column vector and retrieve its 
   * size, i.e. the number of parameters
   */
  if(mxGetM(prhs[0])==1){
    m=mxGetN(prhs[0]);
    isrow_p0 = 1;
  }
  else{
    m=mxGetM(prhs[0]);
    isrow_p0 = 0;
  }
  /* copy input parameter vector to avoid destroying it */
  p=mxMalloc(m*sizeof(double));
  for(i=0; i<m; ++i)
    p[i]=p0[i];
    
  /** x **/
  /* the third required argument must be a real row or column vector */
  if(!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || !(mxGetM(prhs[1])==1 || mxGetN(prhs[1])==1))
    mexErrMsgTxt("fit_front: x must be a real vector.");
  x=mxGetPr(prhs[1]);
  n=__MAX__(mxGetM(prhs[1]), mxGetN(prhs[1]));

  /** itmax **/
  /* the fourth required argument must be a scalar */
  if(!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)
    mexErrMsgTxt("fit_front: itmax must be a scalar.");
  itmax=(int)mxGetScalar(prhs[2]);

  /** stop_tol **/
  /* the fourth required argument must be a scalar */
  if(!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxGetM(prhs[3])!=1 || mxGetN(prhs[3])!=1)
    mexErrMsgTxt("fit_front: itmax must be a scalar.");
  stop_tol=mxGetScalar(prhs[3]);

  opts[1] = stop_tol;
  opts[2] = stop_tol;
  opts[3] = stop_tol;

  /** mintype (optional) **/
  mintype=MIN_CONSTRAINED_BC;

  /* arguments below this point are optional and their presence depends
   * upon the minimization type determined above
   */
  /** lb, ub **/
    /* check if the next two arguments are real row or column vectors */
  if(mxIsDouble(prhs[4]) && !mxIsComplex(prhs[4]) && (mxGetM(prhs[4])==1 || mxGetN(prhs[4])==1)){
    if(mxIsDouble(prhs[5]) && !mxIsComplex(prhs[5]) && (mxGetM(prhs[5])==1 || mxGetN(prhs[5])==1)){
      if((__MAX__(mxGetM(prhs[4]), mxGetN(prhs[4])))!=m)
        mexErrMsgTxt("fit_front: lb must have as many elements as p0");
      if((__MAX__(mxGetM(prhs[5]), mxGetN(prhs[5])))!=m)
        mexErrMsgTxt("fit_front: ub must have as many elements as p0");

      lb=mxGetPr(prhs[4]);
      ub=mxGetPr(prhs[5]);
    }
  }

  /** pos **/
  /* the third required argument must be a real row or column vector */
  if(!mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) || !(mxGetM(prhs[6])==1 || mxGetN(prhs[6])==1) || (__MAX__(mxGetM(prhs[6]), mxGetN(prhs[6])) != n)) {
    mexErrMsgTxt("fit_front: pos must be a real vector.");
  }

  pos=mxGetPr(prhs[6]);

  status=dlevmar_bc_dif(func, p, x, m, n, lb, ub, NULL, itmax, opts, NULL, NULL, NULL, pos);

  /* copy back return results */
  /** ret **/
  plhs[0]=mxCreateDoubleMatrix(1, 1, mxREAL);
  ret=mxGetPr(plhs[0]);
  ret[0]=(double)status;

  /** popt **/
  plhs[1]=(isrow_p0==1)? mxCreateDoubleMatrix(1, m, mxREAL) : mxCreateDoubleMatrix(m, 1, mxREAL);
  pdbl=mxGetPr(plhs[1]);
  for(i=0; i<m; ++i)
    pdbl[i]=p[i];

cleanup:
  /* cleanup */
  mxFree(p);

  if(status==LM_ERROR)
    mexWarnMsgTxt("fit_front: optimization returned with an error!");
}
