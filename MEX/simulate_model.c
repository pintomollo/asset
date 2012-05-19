#include <math.h>
#include <string.h> 
#include "mex.h"

#define MOD(x, y) (int)((x) - (y) * floor((double)(x) / (double)(y)))

static void finite_difference(double *x, double *dx, double h, int npop, int npts, int order) {

  int i, p;
  double *xi, *dxi;
  double x_2, x_1, x0, x1, x2;
  double dnorm, ddnorm;

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

      x_2 = x_1;
      x_1 = x0;
      x0 = x1;
      x1 = x2;
      x2 = xi[npts - abs((i+3) - npts)];
    }
  }

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
    indf = floor(index);
    yf = index - indf;

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

  double *xi;
  int i, p;

  for (p=0; p<npops; p++) {
    xi = x + npts*p;
    ints[p] = 0;
    for (i=1;i<npts-1;i++) {
      ints[p] += 2*xi[i];
    }
    ints[p] += xi[0];
    ints[p] += xi[npts-1];

    ints[p] = ints[p] * (h/2);
  }

  return;
}

/* x0, dt, tmax, */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  double *time, *concentr, *flow, *current_flow, *x, *dx, *ddx, tmax, output_rate, dt, current_dt, t = 0, t_flow, x_step, *x0, *params, *cyto, *x_prev, *tmp, *simulation;
  int ntimes, count = 1, npts, npops, page_size = 500, nflow, nparams, i;
  mwSize m, n;

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

    plhs[0] = mxCreateDoubleMatrix(npops*npts, ntimes, mxREAL);
    simulation = mxGetPr(plhs[0]);

    for (i = 0; i<npts*npops; i++) {
      x[i] = x0[i];
      simulation[i] = x0[i];
    }

    while (t < tmax) {
      bilinear_flow(flow, x, current_flow, t / t_flow, npts, nflow);
      finite_difference(x, ddx, x_step, npops, npts, 2);
      finite_difference(current_flow, dx, x_step, npops, npts, 1);
      trapezoidal_integral(x, x_step, cyto, npts, npops);

      tmp = x_prev;
      x_prev = x;
      x = tmp;
      
      cyto[0] = params[6] - cyto[0] * (params[7] / params[8]);
      cyto[1] = params[6 + nparams] - cyto[1] * (params[7 + nparams] / params[8 + nparams]);

      for (i = 0; i<npts; i++) {
        x[i] = x_prev[i] + dt*(
               params[0]*ddx[i] + 
               params[1]*cyto[0] - 
               params[2]*x_prev[i] - 
               params[3]*pow(x_prev[i], params[4])*pow(x_prev[i + npts], params[5]) - 
               current_flow[i]);

        x[i+npts] = x_prev[i+npts] + dt*(
               params[0 + nparams]*ddx[i+npts] + 
               params[1 + nparams]*cyto[1] - 
               params[2 + nparams]*x_prev[i+npts] - 
               params[3 + nparams]*pow(x_prev[i], params[4 + nparams])*pow(x_prev[i + npts], params[5 + nparams]) - 
               current_flow[i+npts]);
      }

      t += dt;

      if (MOD(t-dt, output_rate) >= MOD(t, output_rate)) {
        memcpy(simulation+count*(npops*npts), x, npops*npts*sizeof(double));
        count++;
      }
    }

    mxFree(x);
    mxFree(dx);
    mxFree(ddx);
    mxFree(x_prev);
    mxFree(cyto);
    mxFree(current_flow);
  }

  return;
}
