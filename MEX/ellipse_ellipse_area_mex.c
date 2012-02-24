#include <math.h>
#include "mex.h"
#include "program_constants.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  double ph1, a1, b1, h1, k1, ph2, a2, b2, h2, k2, val;
  int ret_val = 0, i;
  mxArray *tmp_data;

  if (nrhs != 10) {
    mexErrMsgTxt("Not the correct number of input arguments (10 required) !");
  }

  ph1 = mxGetScalar(prhs[0]);
  a1 = mxGetScalar(prhs[1]);
  b1 = mxGetScalar(prhs[2]);
  h1 = mxGetScalar(prhs[3]);
  k1 = mxGetScalar(prhs[4]);
  ph2 = mxGetScalar(prhs[5]);
  a2 = mxGetScalar(prhs[6]);
  b2 = mxGetScalar(prhs[7]);
  h2 = mxGetScalar(prhs[8]);
  k2 = mxGetScalar(prhs[9]);

  val = (double) ellipse_ellipse_overlap(ph1, a1, b1, h1, k1, ph2, a2, b2, h2, k2, &ret_val);

  plhs[0] = mxCreateDoubleScalar(val);
}
