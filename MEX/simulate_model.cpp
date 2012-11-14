#include <math.h>
#include <string.h>
#include <stdio.h>
#include "mex.h"
#include "rkf45.h"

#define POS(x) ((x) > 0 ? 1 : 0)
#define NEG(x) ((x) < 0 ? 1 : 0)

double *flow, *current_flow, *dx, *ddx, *cyto, *params;
double t_flow, x_step;
int npts, npops, nflow, nparams;

static void finite_difference(double *x, double *dx, double *ddx, double *flow, double h, int npop, int npts) {

  int i, p;
  double *xi, *dxi, *ddxi;
  double x_2, x_1, x0, x1, x2, f_1, f0, f1, cpos, cneg;
  double dnorm, ddnorm;

  dnorm = 1 / h;
  ddnorm = 1 / (12*h*h);

  xi = x;
  dxi = dx;
  ddxi = ddx;

  for (p=0; p<npop; ++p) {

    /* Reflective boundary conditions */
    x_2 = xi[1];
    x_1 = xi[0];
    x0 = xi[0];
    x1 = xi[1];
    x2 = xi[2];

    f_1 = flow[0];
    f0 = flow[0];
    f1 = flow[1];

    for (i=0; i<npts-3; ++i) {
      cpos = POS(f0);
      cneg = NEG(f0);

      dxi[i] = (cpos*(x0*f0 - x_1*f_1) + cneg*(x1*f1 - x0*f0)) * dnorm;
      ddxi[i] = (-x_2 + 16*x_1 - 30*x0 + 16*x1 - x2) * ddnorm;

      x_2 = x_1;
      x_1 = x0;
      x0 = x1;
      x1 = x2;

      x2 = xi[i+3];

      f_1 = f0;
      f0 = f1;
      f1 = flow[i+2];
    }

    cpos = POS(f0);
    cneg = NEG(f0);

    dxi[i] = (cpos*(x0*f0 - x_1*f_1) + cneg*(x1*f1 - x0*f0)) * dnorm;
    ddxi[i] = (-x_2 + 16*x_1 - 30*x0 + 16*x1 - x2) * ddnorm;

    x_2 = x_1;
    x_1 = x0;
    x0 = x1;
    x1 = x2;

    x2 = xi[npts-1];

    f_1 = f0;
    f0 = f1;
    f1 = flow[npts-1];

    i++;

    cpos = POS(f0);
    cneg = NEG(f0);

    dxi[i] = (cpos*(x0*f0 - x_1*f_1) + cneg*(x1*f1 - x0*f0)) * dnorm;
    ddxi[i] = (-x_2 + 16*x_1 - 30*x0 + 16*x1 - x2) * ddnorm;

    x_2 = x_1;
    x_1 = x0;
    x0 = x1;
    x1 = x2;

    x2 = xi[npts-2];

    f_1 = f0;
    f0 = f1;
    f1 = flow[npts-1];

    i++;

    cpos = POS(f0);
    cneg = NEG(f0);

    dxi[i] = (cpos*(x0*f0 - x_1*f_1) + cneg*(x1*f1 - x0*f0)) * dnorm;
    ddxi[i] = (-x_2 + 16*x_1 - 30*x0 + 16*x1 - x2) * ddnorm;

    xi += npts;
    dxi += npts;
    ddxi += npts;
  }

  return;
}

static void bilinear_flow(double *flow, double *current, double index, int npts, int ntimes) {

  double yf, yc;
  int i, indf, indc;

  if (index >= ntimes-1) {
    indf = ntimes-1;
    yf = 0;
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

  for (i=0; i<npts; ++i) {
    current[i] = (flow[i+indf]*yf + flow[i+indc]*yc);
  }

  return;
}

static void trapezoidal_integral(double *x, double h, double *ints, int npts, int npops) {

  double *xi, tmp_int, tmp_x;
  int i, p;

  for (p=0; p<npops; ++p) {
    xi = x + npts*p;
    tmp_int = 0;
    for (i=1;i<npts-1;++i) {
      tmp_x = xi[i];
      tmp_int += 2*tmp_x;
    }

    tmp_int += xi[0];
    tmp_int += xi[npts-1];

    ints[p] = tmp_int * (h/2);
  }

  return;
}

static void reactions(double *fx, double *x_prev, double *dx, double *ddx, double *cyto, double *params, int npts, int npops, int nparams) {

  double *fxi, *xi_prev, *dxi, *ddxi, *paramsi;
  int i;

  fxi = fx;
  xi_prev = x_prev;
  dxi = dx;
  ddxi = ddx;
  paramsi = params;
  
  const double p0 = paramsi[0];
  const double c0 = paramsi[1]*(paramsi[5] - cyto[0] * (paramsi[6] / paramsi[7]));
  const double p2 = paramsi[2];
  const double p3 = paramsi[3];
  const double p4 = paramsi[4];

  for (i=0;i<npts;++i) {
    const double x_p = xi_prev[i];

    fxi[i] = p0*ddxi[i] + 
             c0 - 
             x_p*p2 - 
             x_p*p3*pow(x_prev[i + npts], p4) - 
             dxi[i];
  }

  fxi += npts;
  xi_prev += npts;
  dxi += npts;
  ddxi += npts;
  paramsi += nparams;

  const double p10 = paramsi[0];
  const double c1 = paramsi[1]*(paramsi[5] - cyto[1] * (paramsi[6] / paramsi[7]));
  const double p12 = paramsi[2];
  const double p13 = paramsi[3];
  const double p14 = paramsi[4];

  for (i=0;i<npts;++i) {
    const double x_p = xi_prev[i];

    fxi[i] = p10*ddxi[i] + 
             c1 - 
             x_p*p12 - 
             x_p*p13*pow(x_prev[i], p14) - 
             dxi[i];
  }

  return;
}

static void goehring_step(double t, double *x, double *fx) {

  bilinear_flow(flow, current_flow, t * t_flow, npts, nflow);
  finite_difference(x, dx, ddx, current_flow, x_step, npops, npts);
  trapezoidal_integral(x, x_step, cyto, npts, npops);
  reactions(fx, x, dx, ddx, cyto, params, npts, npops, nparams);
}

/* x0, dt, tmax, */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  double *time, *x, *fx, *x0, *tmp, *simulation, *curr_sim, *timing;
  double tmax, output_rate, dt, current_dt, t = 0, clf_const, inv_x_step, 
        time_count = 0;
  double relerr, abserr;
  int ntimes, count = 1, page_size = 500, i, j, niter, ntotal, ndata, nsteps = 1;
  int flag, prev_flag;
  bool saved = false;
  mwSize m, n;

  if (nrhs != 9) {
    mexErrMsgTxt("Simulation requires 9 input arguments !");
  } else {

    x0 = mxGetPr(prhs[0]);
    npts = mxGetM(prhs[0]);
    npops = mxGetN(prhs[0]);
    ntotal = npts*npops;
    ndata = ntotal*sizeof(double);

    params = mxGetPr(prhs[1]);
    nparams = mxGetM(prhs[1]);
    if (mxGetN(prhs[1]) != npops) {
      mexErrMsgTxt("Parameters have to be defined for every population");
    }

    x_step = mxGetScalar(prhs[2]);
    inv_x_step = 1 / x_step;
    tmax = mxGetScalar(prhs[3]);
    dt = mxGetScalar(prhs[4]);
    output_rate = mxGetScalar(prhs[5]);

    if (dt == 0) {
      mexErrMsgTxt("A time step has to be defined");
    }

    flow = mxGetPr(prhs[6]);
    nflow = mxGetN(prhs[6]);
    if (mxGetM(prhs[6]) != npts) {
      mexErrMsgTxt("Flow has to be defined at every position of the lattice !");
    }

    t_flow = 1 / mxGetScalar(prhs[7]);
    ntimes = ceil(tmax / output_rate) + 1;

    niter = (int)(mxGetScalar(prhs[8]));

    if ((x = (double *)mxCalloc(ntotal, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((fx = (double *)mxCalloc(ntotal*nsteps, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((dx = (double *)mxCalloc(ntotal, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((ddx = (double *)mxCalloc(ntotal, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((current_flow = (double *)mxCalloc(npts, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((cyto = (double *)mxCalloc(npops, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((simulation = (double *)mxCalloc(ntotal*ntimes, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((timing = (double *)mxCalloc(ntimes, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }

    memcpy(simulation, x0, ndata);
    curr_sim = simulation + ntotal;
    memcpy(x, x0, ndata);
    timing[0] = 0;

    abserr = sqrt(r8_epsilon());
    relerr = sqrt(r8_epsilon());

    prev_flag = -1;

    for (j = 0; j < niter; ++j) {
      saved = false;
      current_dt = dt;

      flag = r8_rkf45(goehring_step, ntotal, x, fx, &t, tmax, &relerr, abserr, prev_flag);

      if (flag == 666) {
        mexPrintf("Error due to flag %d\n", prev_flag);
        return;
      } else {
        prev_flag = flag;
      }

      if (t - time_count >= output_rate) {
        saved = true;
        memcpy(curr_sim, x, ndata);
        curr_sim += ntotal;
        timing[count] = t;
        count++;
        time_count += output_rate;

        if (count == ntimes) {
          ntimes += page_size;

          if ((simulation = (double *)mxRealloc(simulation, ntimes * ndata)) == NULL) {
            mexErrMsgTxt("Memory allocation failed !");
          }
          if ((timing = (double *)mxRealloc(timing, ntimes * sizeof(double))) == NULL) {
            mexErrMsgTxt("Memory allocation failed !");
          }
          curr_sim = simulation + (count-1)*ntotal;
        }
      }

      if (t >= tmax) {
        break;
      }
    }
    if (!saved) {
      memcpy(curr_sim, x, ndata);
      timing[count] = t;
    } else {
      count--;
    }
    if (j == niter) {
      mexWarnMsgTxt("Simulation reached the iteration limit !");
    }

    plhs[0] = mxCreateDoubleMatrix(npops*npts, count, mxREAL);
    tmp = mxGetPr(plhs[0]);
    memcpy(tmp, simulation, count*ndata);
    plhs[1] = mxCreateDoubleMatrix(1, count, mxREAL);
    tmp = mxGetPr(plhs[1]);
    memcpy(tmp, timing, count*sizeof(double));

    mxFree(x);
    mxFree(dx);
    mxFree(ddx);
    mxFree(fx);
    mxFree(cyto);
    mxFree(current_flow);
    mxFree(simulation);
    mxFree(timing);
  }

  return;
}
