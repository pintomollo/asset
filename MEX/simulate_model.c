#include <math.h>
#include <string.h>
#include <stdio.h>
#include "mex.h"

#define MOD(x, y) ((x) - (y) * floor((double)(x) / (double)(y)))
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define POS(x) ((x) > 0 ? 1 : 0)
#define NEG(x) ((x) < 0 ? 1 : 0)

static double finite_difference(double *x, double *dx, double *flow, double h, int npop, int npts, int order) {

  int i, p;
  double *xi, *dxi;
  double x_2, x_1, x0, x1, x2, f_1, f0, f1, cpos, cneg;
  double dnorm, ddnorm, max_val = 0;

  dnorm = 1 / h;
  //dnorm = 1 / (12*h);
  ddnorm = 1 / (12*h*h);

  /*for (p=0; p<npop; ++p) {
    xi = x + npts*p;
    dxi = dx + npts*p;

    // Reflective boundary conditions
    x_2 = xi[1];
    x_1 = xi[0];
    x0 = xi[0];
    x1 = xi[1];
    x2 = xi[2];

    for (i=0; i<npts; ++i) {
      if (order == 1) {
        dxi[i] = (x_2 - 8*x_1 + 8*x1 - x2) * dnorm;
      } else {
        dxi[i] = (-x_2 + 16*x_1 - 30*x0 + 16*x1 - x2) * ddnorm;
      }

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
  }*/

  xi = x;
  dxi = dx;

  if (order == 1) {
    for (p=0; p<npop; ++p) {

      /* Reflective boundary conditions */
      //x_2 = xi[1];
      x_1 = xi[0];
      x0 = xi[0];
      x1 = xi[1];

      f_1 = flow[0];
      f0 = flow[0];
      f1 = flow[1];
      //x2 = xi[2];

      for (i=0; i<npts-2; ++i) {
        cpos = POS(f0);
        cneg = NEG(f0);

        dxi[i] = (cpos*(x0*f0 - x_1*f_1) + cneg*(x1*f1 - x0*f0)) * dnorm;
        /*dxi[i] = (x0 - x_1) * dnorm;
        dxi[i] = (x1 - x0) * dnorm;*/
        //dxi[i] = (x_2 - 8*x_1 + 8*x1 - x2) * dnorm;

        if (max_val < fabs(dxi[i])) {
          max_val = fabs(dxi[i]);
        }

        //x_2 = x_1;
        x_1 = x0;
        x0 = x1;
        //x1 = x2;

        x1 = xi[i+2];
        //x2 = xi[i+3];

        f_1 = f0;
        f0 = f1;
        f1 = flow[i+2];
      }

      for (i=npts-2; i<npts; ++i) {
        cpos = POS(f0);
        cneg = NEG(f0);

        dxi[i] = (cpos*(x0*f0 - x_1*f_1) + cneg*(x1*f1 - x0*f0)) * dnorm;
        //dxi[i] = (x_2 - 8*x_1 + 8*x1 - x2) * dnorm;
        //dxi[i] = (x0 - x_1) * dnorm;
        //dxi[i] = (x1 - x0) * dnorm;

        if (max_val < fabs(dxi[i])) {
          max_val = fabs(dxi[i]);
        }

        //x_2 = x_1;
        x_1 = x0;
        x0 = x1;
        //x1 = x2;

        x1 = xi[npts - ((i + 2) - (npts - 1))];
        //x2 = xi[npts - ((i + 3) - (npts - 1))];
        f_1 = f0;
        f0 = f1;
        f1 = flow[npts - ((i + 2) - (npts - 1))];
      }

      xi += npts;
      dxi += npts;
    }
  } else {

    for (p=0; p<npop; ++p) {

      /* Reflective boundary conditions */
      x_2 = xi[1];
      x_1 = xi[0];
      x0 = xi[0];
      x1 = xi[1];
      x2 = xi[2];

      for (i=0; i<npts-3; ++i) {
        dxi[i] = (-x_2 + 16*x_1 - 30*x0 + 16*x1 - x2) * ddnorm;

        x_2 = x_1;
        x_1 = x0;
        x0 = x1;
        x1 = x2;

        x2 = xi[i+3];
      }

      for (i=npts-3; i<npts; ++i) {
        dxi[i] = (-x_2 + 16*x_1 - 30*x0 + 16*x1 - x2) * ddnorm;

        x_2 = x_1;
        x_1 = x0;
        x0 = x1;
        x1 = x2;

        x2 = xi[npts - ((i + 3) - (npts - 1))];
      }

      xi += npts;
      dxi += npts;
    }
  }

  return;
}

static void bilinear_flow(double *flow, double *conc, double *current, double index, int npts, int ntimes) {

  double yf, yc, value;
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
    /*value = (flow[i+indf]*yf + flow[i+indc]*yc);
    current[i] = value * conc[i];
    current[i + npts] = value * conc[i + npts];*/
  }

  return;
}

static void trapezoidal_integral(double *x, double h, double *ints, int npts, int npops) {

  double *xi, count;
  int i, p;

  for (p=0; p<npops; ++p) {
    xi = x + npts*p;
    ints[p] = 0;
    count = 0;
    for (i=1;i<npts-1;++i) {
      ints[p] += 2*xi[i];
      count += xi[i];
    }

    ints[p] += xi[0];
    ints[p] += xi[npts-1];

    ints[p] = ints[p] * (h/2);
  }

  return;
}

static void reactions(double *x, double *x_prev, double *dx, double *ddx, double *cyto, double *params, double dt, int npts, int npops, int nparams) {

  double *xi, *xi_prev, *dxi, *ddxi, *paramsi;
  int i;
  //int i, p, ntotal, incr;

  //ntotal = npts*npops;

  /*for (p=0; p<npops; ++p) {
    incr = npts*p;
    xi = x + incr;
    xi_prev = x_prev + incr;
    dxi = dx + incr;
    ddxi = ddx + incr;
    paramsi = params + nparams*p;

    incr += npts;

    cyto[p] = paramsi[5] - cyto[p] * (paramsi[6] / paramsi[7]);

    for (i=0;i<npts;++i) {
      xi[i] = xi_prev[i] + dt*(
             paramsi[0]*ddxi[i] + 
             paramsi[1]*cyto[p] - 
             paramsi[2]*xi_prev[i] - 
             paramsi[3]*xi_prev[i]*pow(x_prev[(int)MOD(i + incr, ntotal)], paramsi[4]) - 
             dxi[i]);
    }
  }*/

  xi = x;
  xi_prev = x_prev;
  dxi = dx;
  ddxi = ddx;
  paramsi = params;

  //incr += npts;

  cyto[0] = paramsi[5] - cyto[0] * (paramsi[6] / paramsi[7]);

  for (i=0;i<npts;++i) {
    xi[i] = xi_prev[i] + dt*(
           paramsi[0]*ddxi[i] + 
           paramsi[1]*cyto[0] - 
           paramsi[2]*xi_prev[i] - 
           paramsi[3]*xi_prev[i]*pow(x_prev[i + npts], paramsi[4]) - 
           dxi[i]);
  }

  xi += npts;
  xi_prev += npts;
  dxi += npts;
  ddxi += npts;
  paramsi += nparams;

  //incr += npts;

  cyto[1] = paramsi[5] - cyto[1] * (paramsi[6] / paramsi[7]);

  for (i=0;i<npts;++i) {
    xi[i] = xi_prev[i] + dt*(
           paramsi[0]*ddxi[i] + 
           paramsi[1]*cyto[1] - 
           paramsi[2]*xi_prev[i] - 
           paramsi[3]*xi_prev[i]*pow(x_prev[i], paramsi[4]) - 
           dxi[i]);
  }

  return;
}

/* x0, dt, tmax, */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  double *time, *concentr, *flow, *current_flow, *x, *dx, *ddx, tmax, output_rate,
         dt, current_dt, t = 0, t_flow, x_step, *x0, *params, *cyto, *x_prev, *tmp,
         *simulation, *curr_sim, *timing, clf_const, inv_x_step, time_count = 0;
  int ntimes, count = 1, npts, npops, page_size = 500, nflow, nparams, i, j, niter, ntotal, ndata;
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

    if (MAX(params[0], params[0+npts])*dt*3 / pow(x_step, 2) >= 1) {
      dt = pow(x_step, 2) / (3*MAX(params[0], params[0+nparams]));
    }

    flow = mxGetPr(prhs[6]);
    nflow = mxGetN(prhs[6]);
    if (mxGetM(prhs[6]) != npts) {
      mexErrMsgTxt("Flow has to be defined at every position of the lattice !");
    }

    t_flow = 1 / mxGetScalar(prhs[7]);
    ntimes = ceil(tmax / output_rate) + 1;

    niter = (int)(mxGetScalar(prhs[8]));

    if ((x = mxCalloc(ntotal, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((x_prev = mxCalloc(ntotal, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((dx = mxCalloc(ntotal, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((ddx = mxCalloc(ntotal, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((current_flow = mxCalloc(npts, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    /*if ((current_flow = mxCalloc(ntotal, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }*/
    if ((cyto = mxCalloc(npops, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((simulation = mxCalloc(ntotal*ntimes, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((timing = mxCalloc(ntimes, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }

    memcpy(simulation, x0, ndata);
    curr_sim = simulation + ntotal;
    memcpy(x, x0, ndata);
    timing[0] = 0;

    for (j = 0; j < niter; ++j) {
      saved = false;
      current_dt = dt;

      bilinear_flow(flow, x, current_flow, t * t_flow, npts, nflow);
      clf_const = finite_difference(x, dx, current_flow, x_step, npops, npts, 1);
      finite_difference(x, ddx, current_flow, x_step, npops, npts, 2);
      trapezoidal_integral(x, x_step, cyto, npts, npops);

      if (clf_const*dt * inv_x_step >= 0.5) {
        current_dt = 0.5 * x_step / clf_const;
      }

      tmp = x_prev;
      x_prev = x;
      x = tmp;

      reactions(x, x_prev, dx, ddx, cyto, params, current_dt, npts, npops, nparams);

      t += current_dt;
      time_count += current_dt;

      if (time_count >= output_rate) {
        saved = true;
        memcpy(curr_sim, x, ndata);
        curr_sim += ntotal;
        timing[count] = t;
        count++;
        time_count -= output_rate;

        if (count == ntimes) {
          ntimes += page_size;

          if ((simulation = mxRealloc(simulation, ntimes * ndata)) == NULL) {
            mexErrMsgTxt("Memory allocation failed !");
          }
          if ((timing = mxRealloc(timing, ntimes * sizeof(double))) == NULL) {
            mexErrMsgTxt("Memory allocation failed !");
          }
          curr_sim = simulation + (count-1)*ntotal;
        }
      }

      if (t > tmax) {
        break;
      }
    }
    if (~saved) {
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
    mxFree(x_prev);
    mxFree(cyto);
    mxFree(current_flow);
    mxFree(simulation);
    mxFree(timing);
  }

  return;
}
