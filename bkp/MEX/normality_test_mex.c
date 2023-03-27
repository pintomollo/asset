#include <math.h>
#include "mex.h"

#ifdef _MSC_VER
double erf(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    double t, y;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    t = 1.0/(1.0 + p*x);
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}
#endif

#define MIN(a,b) ((a) > (b) ? (b) : (a))
#define MAX(a,b) ((a) < (b) ? (b) : (a))
#define POW2(a) ((a)*(a))

/* Compute the normality test using the Kolmogorov-Smirnov test */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  double *x, *y, *mu, *cdf, *bounds, *ret;
  double center = 0.0, x_min, x_mu, x_max, y_min, y_max, x_curr, x_prev, width;
  double c_min, c_max, x_cdf, y_cdf, slope, bkg, score, best_score, overall_best;
  double best_peak;
  int *bound_indexes, *widths;
  int n, m, indx, w_min, w_max, i_min, i_max, coeff = 3, max_width = 0, best_width;
  int overall_width;
  register int i, j;

  if (nrhs < 2) {
    mexErrMsgTxt("Not enough input arguments (2 is the minimum) !");
  } 
  
  x = mxGetPr(prhs[0]);
  y = mxGetPr(prhs[1]);

  n = mxGetNumberOfElements(prhs[0]);
  if (mxGetNumberOfElements(prhs[1]) != n) {
    mexErrMsgTxt("x and y should be vectors of the same size !");
  }

  if (nrhs > 2) {

    mu = mxGetPr(prhs[2]);
    m = mxGetNumberOfElements(prhs[2]);

    if (nrhs > 3) {
      coeff = (int)mxGetScalar(prhs[3]);
    }
  } else {
    mu = &center;
    m = 1;
  }

  if ((bound_indexes = (int*)mxCalloc(3*m, sizeof(int))) == NULL) {
    mexErrMsgTxt("Memory allocation failed !");
  }

  if ((bounds = (double*)mxMalloc(2*m*sizeof(double))) == NULL) {
    mexErrMsgTxt("Memory allocation failed !");
  }

  if ((widths = (int*)mxMalloc(2*m*sizeof(int))) == NULL) {
    mexErrMsgTxt("Memory allocation failed !");
  }

  if ((cdf = (double*)mxMalloc(n*sizeof(double))) == NULL) {
    mexErrMsgTxt("Memory allocation failed !");
  }

  /* Compute the possible sizes of gaussians */
  x_min = x[1];
  x_max = x[n-2];
  for (i=0; i < m; i++) {
    x_curr = mu[i];
    width = floor(MIN(x_curr - x_min, x_max - x_curr)/coeff);
    widths[2*i] = width;
    width = width*coeff;
    bounds[2*i] = x_curr - width;
    bounds[2*i+1] = x_curr + width;
  }

  /* The CDF of the input signal */
  cdf[0] = 0;
  for (i=1; i < n; i++) {
    x_curr = x[i];
    x_prev = x[i-1];
    cdf[i] = cdf[i-1] + 0.5*(y[i]+y[i-1]);

    for (j=0; j < m; j++) {
      if (x_curr == bounds[2*j]) {
        bound_indexes[3*j] = i-1;
      } else if (x_curr == mu[j]) {
        bound_indexes[3*j+1] = i;
      } else if (x_curr == bounds[2*j+1]) {
        bound_indexes[3*j+2] = i+1;
      }
    }
  }

  /* The last bounds */
  for (i=0; i < m; i++) {
    indx = bound_indexes[3*i+1];

    if (indx < 3 || indx > (n-4)) {
      widths[2*i+1] = widths[2*i]+1;
    } else {
      widths[2*i+1] = 2;
    }
  }

  overall_best = 1;
  overall_width = -1;
  best_peak = 0;

  for (i=0; i < m; i++) {
    w_max = widths[2*i];
    w_min = widths[2*i+1];
    best_score = 1;
    best_width = -1;
    for (j=w_min; j<=w_max; j++) {
      width = coeff*j;
      x_mu = mu[i];
      i_min = bound_indexes[3*i+1]-width;
      i_max = bound_indexes[3*i+1]+width;
      c_min = cdf[i_min-1];
      c_max = cdf[i_max+1];
      x_min = x[i_min-1];
      bkg = y[i_min-1];
      slope = 0.25*(y[i_max+1]-bkg) / (width+1);
      c_max = c_max - c_min - (slope*POW2(2*(width+1)) + bkg*(2*(width+1)));
      c_max = 1/c_max;
      width = 1/(sqrt(2)*(width/coeff));

      score = 0;
      for (n=i_min; n<=i_max; n++) {
        x_curr = x[n];
        x_cdf = 0.5 + 0.5*erf((x[n] - x_mu) * width);
        x_curr -= x_min;
        y_cdf = (cdf[n] - c_min - (slope*POW2(x_curr) + bkg*(x_curr)))*c_max;
        score = MAX(score, fabs(y_cdf-x_cdf));
      }

      if (best_score > score) {
        best_score = score;
        best_width = j;
      }
    }
    if (overall_best > best_score) {
      overall_best = best_score;
      overall_width = best_width;
      best_peak = mu[i];
    }
  }

  plhs[0]=mxCreateDoubleMatrix(1, 2, mxREAL);
  ret=mxGetPr(plhs[0]);
  ret[0] = best_peak;
  ret[1] = (double)overall_width;

  mxFree(bounds);
  mxFree(widths);
  mxFree(bound_indexes);
  mxFree(cdf);

  return;
}
