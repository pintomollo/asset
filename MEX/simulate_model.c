#include <math.h>
#include <string.h>
#include <stdio.h>
#include "mex.h"

#define MOD(x, y) ((x) - (y) * floor((double)(x) / (double)(y)))
#define max(x, y) ((x) > (y) ? (x) : (y))

static double finite_difference(double *x, double *dx, double h, int npop, int npts, int order) {

  int i, p;
  double *xi, *dxi;
  double x_2, x_1, x0, x1, x2;
  double dnorm, ddnorm, max_val = 0;

  dnorm = 12*h;
  ddnorm = dnorm*h;

  for (p=0; p<npop; p++) {
    xi = x + npts*p;
    dxi = dx + npts*p;

    /* Reflective boundary conditions */
    x_2 = xi[1];
    x_1 = xi[0];
    x0 = xi[0];
    x1 = xi[1];
    x2 = xi[2];

    for (i=0; i<npts; i++) {
      if (order == 1) {
        dxi[i] = (x_2 - 8*x_1 + 8*x1 - x2) / dnorm;
      } else {
        dxi[i] = (-x_2 + 16*x_1 - 30*x0 + 16*x1 - x2) / ddnorm;
      }

      //if (p==1)
      //  fprintf(fid, "%f %f\n", dx[i], dxi[i]);

      if (max_val < fabs(dxi[i])) {
        max_val = fabs(dxi[i]);
      }

      x_2 = x_1;
      x_1 = x0;
      x0 = x1;
      x1 = x2;

      if (i + 3 >= npts) {
        x2 = xi[npts - ((i + 3) - (npts - 1))];
      } else {
        x2 = xi[i+3];
      }
    }
  }
  //fprintf(fid, "........................\n");

  return;
}

static void bilinear_flow(double *flow, double *conc, double *current, double index, int npts, int ntimes) {

  double yf, yc, value;
  int i, indf, indc;

  if (index >= ntimes-1) {
    indf = ntimes-1;
    yf = 1;
    indc = ntimes-1;
    yc = 0;
  } else {
    indf = (int)floor(index);
    yf = 1 - (index - (double)indf);

    if (yf == 0) {
      indc = indf;
      yc = 1;
    } else {
      indc = indf + 1;
      yc = 1 - yf;
    }
  }

  indf = indf*npts;
  indc = indc*npts;

  for (i=0; i<npts; i++) {
    value = (flow[i+indf]*yf + flow[i+indc]*yc);
    current[i] = value * conc[i];
    current[i + npts] = value * conc[i + npts];
  }

  return;
}

static void trapezoidal_integral(double *x, double h, double *ints, int npts, int npops) {

  double *xi, count;
  int i, p;

  for (p=0; p<npops; p++) {
    xi = x + npts*p;
    ints[p] = 0;
    count = 0;
    for (i=1;i<npts-1;i++) {
      ints[p] += 2*xi[i];
      count += xi[i];
    }

    ints[p] += xi[0];
    ints[p] += xi[npts-1];

    ints[p] = ints[p] * (h/2);
  }

  return;
}

/* x0, dt, tmax, */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  double *time, *concentr, *flow, *current_flow, *x, *dx, *ddx, tmax, output_rate,
         dt, current_dt, t = 0, t_flow, x_step, *x0, *params, *cyto, *x_prev, *tmp,
         *simulation, *timing, clf_const;
  int ntimes, count = 1, npts, npops, page_size = 500, nflow, nparams, i;
  bool saved = false;
  mwSize m, n;
  FILE *fid;

  if (nrhs != 8) {
    mexErrMsgTxt("Simulation requires 8 input arguments !");
  } else {

    x0 = mxGetPr(prhs[0]);
    npts = mxGetM(prhs[0]);
    npops = mxGetN(prhs[0]);

    params = mxGetPr(prhs[1]);
    nparams = mxGetM(prhs[1]);
    if (mxGetN(prhs[1]) != npops) {
      mexErrMsgTxt("Parameters have to be defined for every population");
    }

    x_step = mxGetScalar(prhs[2]);
    tmax = mxGetScalar(prhs[3]);
    dt = mxGetScalar(prhs[4]);
    output_rate = mxGetScalar(prhs[5]);

    if (max(params[1], params[1+npts])*dt*3 / pow(x_step, 2) >= 1) {
      dt = pow(x_step, 2) / (3*max(params[1], params[1+npts]));
    }

    flow = mxGetPr(prhs[6]);
    nflow = mxGetN(prhs[6]);
    if (mxGetM(prhs[6]) != npts) {
      mexErrMsgTxt("Flow has to be defined at every position of the lattice !");
    }

    t_flow = mxGetScalar(prhs[7]);
    ntimes = ceil(tmax / output_rate) + 1;

    if ((x = mxCalloc(npts*npops, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((x_prev = mxCalloc(npts*npops, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((dx = mxCalloc(npts*npops, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((ddx = mxCalloc(npts*npops, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((current_flow = mxCalloc(npts*npops, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((cyto = mxCalloc(npops, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((simulation = mxCalloc(npops*npts*ntimes, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((timing = mxCalloc(ntimes, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }

    memcpy(simulation, x0, npops*npts*sizeof(double));
    memcpy(x, x0, npops*npts*sizeof(double));
    timing[0] = 0;

  //fid = fopen("data1.txt", "w");

    while (t < tmax) {

//      mexPrintf("%f\n", t);

      saved = false;
      current_dt = dt;

      bilinear_flow(flow, x, current_flow, t / t_flow, npts, nflow);
      clf_const = finite_difference(current_flow, dx, x_step, npops, npts, 1);
      finite_difference(x, ddx, x_step, npops, npts, 2);
      trapezoidal_integral(x, x_step, cyto, npts, npops);

      if (clf_const*dt / x_step >= 0.5) {
        current_dt = 0.5 * x_step / clf_const;
      }

      tmp = x_prev;
      x_prev = x;
      x = tmp;

      cyto[0] = params[5] - cyto[0] * (params[6] / params[7]);
      cyto[1] = params[5 + nparams] - cyto[1] * (params[6 + nparams] / params[7 + nparams]);

      for (i = 0; i<npts; i++) {
        x[i] = x_prev[i] + current_dt*(
               params[0]*ddx[i] + 
               params[1]*cyto[0] - 
               params[2]*x_prev[i] - 
               params[3]*x_prev[i]*pow(x_prev[i + npts], params[4]) - 
               dx[i]);

        x[i+npts] = x_prev[i+npts] + current_dt*(
               params[0 + nparams]*ddx[i+npts] + 
               params[1 + nparams]*cyto[1] - 
               params[2 + nparams]*x_prev[i+npts] - 
               params[3 + nparams]*pow(x_prev[i], params[4 + nparams])*x_prev[i + npts] - 
               dx[i+npts]);

/*        fprintf(fid, "%f %f\n", params[1]*cyto[0] - 
               params[2]*x_prev[i] - 
               params[3]*x_prev[i]*pow(x_prev[i + npts], params[4]),
               params[1 + nparams]*cyto[1] - 
               params[2 + nparams]*x_prev[i+npts] - 
               params[3 + nparams]*pow(x_prev[i], params[4 + nparams])*x_prev[i + npts]);*/

      }
 //     fprintf(fid, "||||||||||||||||||||||||\n");

      t += current_dt;

      if (MOD(t-current_dt, output_rate) >= MOD(t, output_rate)) {
        saved = true;
        memcpy(simulation+count*(npops*npts), x, npops*npts*sizeof(double));
        timing[count] = t;
        count++;

        if ((count % ntimes) == 0) {
          ntimes += page_size;

          if ((simulation = mxRealloc(simulation, ntimes * npops*npts*sizeof(double))) == NULL) {
            mexErrMsgTxt("Memory allocation failed !");
          }
          if ((timing = mxRealloc(timing, ntimes * sizeof(double))) == NULL) {
            mexErrMsgTxt("Memory allocation failed !");
          }
        }
      }
    }
    if (~saved) {
      memcpy(simulation+count*(npops*npts), x, npops*npts*sizeof(double));
      timing[count] = t;
    } else {
      count--;
    }

    plhs[0] = mxCreateDoubleMatrix(npops*npts, count, mxREAL);
    tmp = mxGetPr(plhs[0]);
    memcpy(tmp, simulation, npops*npts*count*sizeof(double));
    plhs[1] = mxCreateDoubleMatrix(1, count, mxREAL);
    tmp = mxGetPr(plhs[1]);
    memcpy(tmp, timing, count*sizeof(double));

    mxFree(x);
    mxFree(dx);
    mxFree(ddx);
    mxFree(x_prev);
    mxFree(cyto);
    mxFree(current_flow);
    mxFree(simulation);
    mxFree(timing);

    //fclose(fid);
  }

  return;
}
