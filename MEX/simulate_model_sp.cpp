#include <xmmintrin.h> //SSE
#include <emmintrin.h> //SSE2
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "mex.h"
#include "rkf45.h"

#define POS(x) ((x) > 0 ? 1.0f : 0.0f)
#define NEG(x) ((x) < 0 ? 1.0f : 0.0f)

#define FASTPOWD(x, y) exp2d2(_mm_mul_pd(log2d2(x), (y)))
#define FASTPOWS(x, y) exp2f4(_mm_mul_ps(log2f4(x), (y)))
#define FASTPOWDS(x, y) _mm_cvtps_pd(exp2f4(_mm_mul_ps(log2f4(_mm_cvtpd_ps(x)), _mm_cvtpd_ps(y))))

#define EXP_POLY_DEGREE 3
#define LOG_POLY_DEGREE 5

#define POLY0(x, c0) _mm_set1_ps(c0)
#define POLY1(x, c0, c1) _mm_add_ps(_mm_mul_ps(POLY0(x, c1), x), _mm_set1_ps(c0))
#define POLY2(x, c0, c1, c2) _mm_add_ps(_mm_mul_ps(POLY1(x, c1, c2), x), _mm_set1_ps(c0))
#define POLY3(x, c0, c1, c2, c3) _mm_add_ps(_mm_mul_ps(POLY2(x, c1, c2, c3), x), _mm_set1_ps(c0))
#define POLY4(x, c0, c1, c2, c3, c4) _mm_add_ps(_mm_mul_ps(POLY3(x, c1, c2, c3, c4), x), _mm_set1_ps(c0))
#define POLY5(x, c0, c1, c2, c3, c4, c5) _mm_add_ps(_mm_mul_ps(POLY4(x, c1, c2, c3, c4, c5), x), _mm_set1_ps(c0))

#define POLYD0(x, c0) _mm_set1_pd(c0)
#define POLYD1(x, c0, c1) _mm_add_pd(_mm_mul_pd(POLYD0(x, c1), x), _mm_set1_pd(c0))
#define POLYD2(x, c0, c1, c2) _mm_add_pd(_mm_mul_pd(POLYD1(x, c1, c2), x), _mm_set1_pd(c0))
#define POLYD3(x, c0, c1, c2, c3) _mm_add_pd(_mm_mul_pd(POLYD2(x, c1, c2, c3), x), _mm_set1_pd(c0))
#define POLYD4(x, c0, c1, c2, c3, c4) _mm_add_pd(_mm_mul_pd(POLYD3(x, c1, c2, c3, c4), x), _mm_set1_pd(c0))
#define POLYD5(x, c0, c1, c2, c3, c4, c5) _mm_add_pd(_mm_mul_pd(POLYD4(x, c1, c2, c3, c4, c5), x), _mm_set1_pd(c0))

inline __m128 exp2f4(__m128 x)
{
  __m128i ipart;
  __m128 fpart, expipart, expfpart;

  x = _mm_min_ps(x, _mm_set1_ps( 129.00000f));
  x = _mm_max_ps(x, _mm_set1_ps(-126.99999f));

  /* ipart = int(x - 0.5) */
  ipart = _mm_cvtps_epi32(_mm_sub_ps(x, _mm_set1_ps(0.5f)));

  /* fpart = x - ipart */
  fpart = _mm_sub_ps(x, _mm_cvtepi32_ps(ipart));

  /* expipart = (float) (1 << ipart) */
  expipart = _mm_castsi128_ps(_mm_slli_epi32(_mm_add_epi32(ipart, _mm_set1_epi32(127)), 23));

  /* minimax polynomial fit of 2**x, in range [-0.5, 0.5[ */
#if EXP_POLY_DEGREE == 5
  expfpart = POLY5(fpart, 9.9999994e-1f, 6.9315308e-1f, 2.4015361e-1f, 5.5826318e-2f, 8.9893397e-3f, 1.8775767e-3f);
#elif EXP_POLY_DEGREE == 4
  expfpart = POLY4(fpart, 1.0000026f, 6.9300383e-1f, 2.4144275e-1f, 5.2011464e-2f, 1.3534167e-2f);
#elif EXP_POLY_DEGREE == 3
  expfpart = POLY3(fpart, 9.9992520e-1f, 6.9583356e-1f, 2.2606716e-1f, 7.8024521e-2f);
#elif EXP_POLY_DEGREE == 2
  expfpart = POLY2(fpart, 1.0017247f, 6.5763628e-1f, 3.3718944e-1f);
#else
#error
#endif

  return _mm_mul_ps(expipart, expfpart);
}

inline __m128 log2f4(__m128 x)
{
  __m128i exp = _mm_set1_epi32(0x7F800000);
  __m128i mant = _mm_set1_epi32(0x007FFFFF);

  __m128 one = _mm_set1_ps( 1.0f);

  __m128i i = _mm_castps_si128(x);

  __m128 e = _mm_cvtepi32_ps(_mm_sub_epi32(_mm_srli_epi32(_mm_and_si128(i, exp), 23), _mm_set1_epi32(127)));

  __m128 m = _mm_or_ps(_mm_castsi128_ps(_mm_and_si128(i, mant)), one);

  __m128 p;

  /* Minimax polynomial fit of log2(x)/(x - 1), for x in range [1, 2[ */
#if LOG_POLY_DEGREE == 6
  p = POLY5( m, 3.1157899f, -3.3241990f, 2.5988452f, -1.2315303f,  3.1821337e-1f, -3.4436006e-2f);
#elif LOG_POLY_DEGREE == 5
  p = POLY4(m, 2.8882704548164776201f, -2.52074962577807006663f, 1.48116647521213171641f, -0.465725644288844778798f, 0.0596515482674574969533f);
#elif LOG_POLY_DEGREE == 4
  p = POLY3(m, 2.61761038894603480148f, -1.75647175389045657003f, 0.688243882994381274313f, -0.107254423828329604454f);
#elif LOG_POLY_DEGREE == 3
  p = POLY2(m, 2.28330284476918490682f, -1.04913055217340124191f, 0.204446009836232697516f);
#else
#error
#endif

  /* This effectively increases the polynomial degree by one, but ensures that log2(1) == 0*/
  p = _mm_mul_ps(p, _mm_sub_ps(m, one));

  return _mm_add_ps(p, e);
}

float *flow, *current_flow, *dx, *ddx, *cyto, *params;
float t_flow, x_step;
int npts, npops, nflow, nparams;
FILE *fid;

static void finite_difference(float *x, float *dx, float *ddx, float *flow, float h, int npop, int npts) {

  int i, p;
  float *xi, *dxi, *ddxi;
  float x_2, x_1, x0, x1, x2, f_2, f_1, f0, f1, f2, cpos, cneg;
  float dnorm, ddnorm;

  dnorm = 1 / h;
  ddnorm = 1 / (12*h*h);

  xi = x;
  dxi = dx;
  ddxi = ddx;

//  fprintf(fid, "%e %e %e %e %e\n", flow[0], flow[1], flow[2], flow[3], flow[4]);
//  fprintf(fid, "%e %e %e %e %e\n", x[0], x[1], x[2], x[3], x[4]);

  for (p=0; p<npop; ++p) {

    /* Reflective boundary conditions */
    x_2 = xi[1];
    x_1 = xi[0];
    x0 = xi[0];
    x1 = xi[1];
    x2 = xi[2];

    f_2 = flow[1];
    f_1 = flow[0];
    f0 = flow[0];
    f1 = flow[1];
    f2 = flow[2];

    for (i=0; i<npts-3; ++i) {
      cpos = POS(f0);
      cneg = NEG(f0);

      dxi[i] = (cpos*((3*(x0*f0) - 4*(x_1*f_1))+(x_2*f_2)) + cneg*(-(3*(x0*f0) + (x2*f2)) + 4*(x1*f1))) * dnorm;
      //ddxi[i] = (-x_2 + 16*x_1 - 30*x0 + 16*x1 - x2) * ddnorm;
      ddxi[i] = ((16*(x_1+x1) - (x2+x_2)) - 30*x0) * ddnorm;

      /*if ((i%2) == 0) {
      fprintf(fid, "%f %f (? %e) (? %e) (? %e) (? %e) (? %e)\n", cpos, cneg, x_2*f_2, x_1*f_1, x0*f0, x1*f1, x2*f2);
      fprintf(fid, "%e %e %e %e %e; %e %e %e %e\n", 4*(x_1*f_1), 3*(x0*f0), (3*(x0*f0)) - (4*(x_1*f_1)), (3*(x0*f0) - 4*(x_1*f_1)) + (x_2*f_2), cpos*((3*(x0*f0) - 4*(x_1*f_1)) + (x_2*f_2)), 3*(x0*f0), 3*(x0*f0) + (x2*f2), (-(3*(x0*f0) +(x2*f2)) + 4*(x1*f1)), cneg*(-(3*(x0*f0) +(x2*f2)) + 4*(x1*f1)));
      //fprintf(fid, "%f (%e %e) %e (%e %e) %e (%e %e) %e : %e\n", cneg, x0, f0, 3*x0*f0, x2, f2, 3*x0*f0 + x2*f2, x1, f1, cneg*(-(3*x0*f0 +x2*f2) - 4*x_1*f_1));
      }*/

      //fprintf(fid, "%e\n%e\n%e\n%e\n%e\n", x_2, x_1, x0, x1, x2);

      x_2 = x_1;
      x_1 = x0;
      x0 = x1;
      x1 = x2;

      x2 = xi[i+3];

      f_2 = f_1;
      f_1 = f0;
      f0 = f1;
      f1 = f2;

      f2 = flow[i+3];
    }

    cpos = POS(f0);
    cneg = NEG(f0);

    //fprintf(fid, "%d)\n", i);
    //dxi[i] = (cpos*(x0*f0 - x_1*f_1) + cneg*(x1*f1 - x0*f0)) * dnorm;
    dxi[i] = (cpos*((3*(x0*f0) - 4*(x_1*f_1))+(x_2*f_2)) + cneg*(-(3*(x0*f0) + (x2*f2)) + 4*(x1*f1))) * dnorm;
    //dxi[i] = (cpos*(3*x0*f0 - 4*x_1*f_1+x_2*f_2) + cneg*(-(3*x0*f0 + x2*f2) + 4*x1*f1)) * dnorm;
    //dxi[i] = (cpos*(3*x0*f0 - 4*x_1*f_1+x_2*f_2) + cneg*(-x2*f2 + 4*x1*f1 - 3*x0*f0)) * dnorm;
    //ddxi[i] = (-x_2 + 16*x_1 - 30*x0 + 16*x1 - x2) * ddnorm;
    ddxi[i] = ((16*(x_1+x1) - (x2+x_2)) - 30*x0) * ddnorm;

    /*  fprintf(fid, "%f %f (? %e) (? %e) (? %e) (? %e) (? %e)\n", cpos, cneg, x_2*f_2, x_1*f_1, x0*f0, x1*f1, x2*f2);
      fprintf(fid, "%e %e %e %e %e; %e %e %e %e\n", 4*(x_1*f_1), 3*(x0*f0), (3*(x0*f0)) - (4*(x_1*f_1)), (3*(x0*f0) - 4*(x_1*f_1)) + (x_2*f_2), cpos*((3*(x0*f0) - 4*(x_1*f_1)) + (x_2*f_2)), 3*(x0*f0), 3*(x0*f0) + (x2*f2), (-(3*(x0*f0) +(x2*f2)) + 4*(x1*f1)), cneg*(-(3*(x0*f0) +(x2*f2)) + 4*(x1*f1)));*/

    x_2 = x_1;
    x_1 = x0;
    x0 = x1;
    x1 = x2;

    x2 = xi[npts-1];

    f_2 = f_1;
    f_1 = f0;
    f0 = f1;
    f1 = f2;
    f2 = flow[npts-1];

    i++;

    cpos = POS(f0);
    cneg = NEG(f0);

//    fprintf(fid, "%d)\n", i);
    //dxi[i] = (cpos*(x0*f0 - x_1*f_1) + cneg*(x1*f1 - x0*f0)) * dnorm;
    //dxi[i] = (cpos*(3*x0*f0 - 4*x_1*f_1+x_2*f_2) + cneg*(-x2*f2 + 4*x1*f1 - 3*x0*f0)) * dnorm;
    //dxi[i] = (cpos*(3*x0*f0 - 4*x_1*f_1+x_2*f_2) + cneg*(-(3*x0*f0 + x2*f2) + 4*x1*f1)) * dnorm;
    dxi[i] = (cpos*((3*(x0*f0) - 4*(x_1*f_1))+(x_2*f_2)) + cneg*(-(3*(x0*f0) + (x2*f2)) + 4*(x1*f1))) * dnorm;
    //ddxi[i] = (-x_2 + 16*x_1 - 30*x0 + 16*x1 - x2) * ddnorm;
    ddxi[i] = ((16*(x_1+x1) - (x2+x_2)) - 30*x0) * ddnorm;

/*      fprintf(fid, "%f %f (? %e) (? %e) (? %e) (? %e) (? %e)\n", cpos, cneg, x_2*f_2, x_1*f_1, x0*f0, x1*f1, x2*f2);
      fprintf(fid, "%e %e %e %e %e; %e %e %e %e\n", 4*(x_1*f_1), 3*(x0*f0), (3*(x0*f0)) - (4*(x_1*f_1)), (3*(x0*f0) - 4*(x_1*f_1)) + (x_2*f_2), cpos*((3*(x0*f0) - 4*(x_1*f_1)) + (x_2*f_2)), 3*(x0*f0), 3*(x0*f0) + (x2*f2), (-(3*(x0*f0) +(x2*f2)) + 4*(x1*f1)), cneg*(-(3*(x0*f0) +(x2*f2)) + 4*(x1*f1)));
    fprintf(fid, "%e ", (cpos*((3*(x0*f0) - 4*(x_1*f_1))+(x_2*f_2)) + cneg*(-(3*(x0*f0) + (x2*f2)) + 4*(x1*f1))));
    fprintf(fid, "%e\n", (cpos*((3*(x0*f0) - 4*(x_1*f_1))+(x_2*f_2)) + cneg*(-(3*(x0*f0) + (x2*f2)) + 4*(x1*f1))) * dnorm);*/

    x_2 = x_1;
    x_1 = x0;
    x0 = x1;
    x1 = x2;

    x2 = xi[npts-2];

    f_2 = f_1;
    f_1 = f0;
    f0 = f1;
    f1 = f2;

    f2 = flow[npts-2];

    i++;

    cpos = POS(f0);
    cneg = NEG(f0);

    //fprintf(fid, "%d)\n", i);
    //dxi[i] = (cpos*(x0*f0 - x_1*f_1) + cneg*(x1*f1 - x0*f0)) * dnorm;
    //dxi[i] = (cpos*(3*x0*f0 - 4*x_1*f_1+x_2*f_2) + cneg*(-(3*x0*f0 + x2*f2) + 4*x1*f1)) * dnorm;
    dxi[i] = (cpos*((3*(x0*f0) - 4*(x_1*f_1))+(x_2*f_2)) + cneg*(-(3*(x0*f0) + (x2*f2)) + 4*(x1*f1))) * dnorm;
    //dxi[i] = (cpos*(3*x0*f0 - 4*x_1*f_1+x_2*f_2) + cneg*(4*x1*f1 - (3*x0*f0 + x2*f2))) * dnorm;
    //ddxi[i] = (-x_2 + 16*x_1 - 30*x0 + 16*x1 - x2) * ddnorm;
    ddxi[i] = ((16*(x_1+x1) - (x2+x_2)) - 30*x0) * ddnorm;

      /*fprintf(fid, "%f %f (? %e) (? %e) (? %e) (? %e) (? %e)\n", cpos, cneg, x_2*f_2, x_1*f_1, x0*f0, x1*f1, x2*f2);
      fprintf(fid, "%e %e %e %e %e; %e %e %e %e\n", 4*(x_1*f_1), 3*(x0*f0), (3*(x0*f0)) - (4*(x_1*f_1)), (3*(x0*f0) - 4*(x_1*f_1)) + (x_2*f_2), cpos*((3*(x0*f0) - 4*(x_1*f_1)) + (x_2*f_2)), 3*(x0*f0), 3*(x0*f0) + (x2*f2), (-(3*(x0*f0) +(x2*f2)) + 4*(x1*f1)), cneg*(-(3*(x0*f0) +(x2*f2)) + 4*(x1*f1)));*/

    xi += npts;
    dxi += npts;
    ddxi += npts;
  }

 /*for (i=0;i<npts;i++){
    fprintf(fid, "%f %f: %.12e %.12e, %f %f: %.12e %.12e\n", x[i], flow[i], dx[i], ddx[i], x[i+npts], flow[i], dx[i+npts], ddx[i+npts]);
  }*/


  return;
}

static void bilinear_flow(float *flow, float *current, float index, int npts, int ntimes) {

  float yf, yc;
  int i, indf, indc;

  if (index >= ntimes-1) {
    indf = ntimes-1;
    yf = 0;
    indc = ntimes-1;
    yc = 0;
  } else {
    indf = (int)floor(index);
    yf = 1 - (index - (float)indf);

    if (yf == 0) {
      indc = indf;
      yc = 1;
    } else {
      indc = indf + 1;
      yc = 1 - yf;
    }
  }
  //fprintf(fid, "%f: %d %f %d %f\n", index, indf, yf, indc, yc);

  indf = indf*npts;
  indc = indc*npts;

  for (i=0; i<npts; ++i) {
    current[i] = (flow[i+indf]*yf + flow[i+indc]*yc);
    /*fprintf(fid, "(%e %e) %e ", flow[i+indf], flow[i+indc], current[i]);

    if ((i%4)==3){
      fprintf(fid, "\n");
    }*/
  }

  return;
}

static void trapezoidal_integral(float *x, float h, float *ints, int npts, int npops) {

  float *xi, tmp_int, tmp_x;
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

static void reactions(float *fx, float *x_prev, float *dx, float *ddx, float *cyto, float *params, int npts, int npops, int nparams) {

  float *paramsi;
  float *tests;

  paramsi = params;

  __m128 p00 = _mm_set1_ps(paramsi[0]);
  __m128 c00 = _mm_set1_ps(paramsi[1]*(paramsi[5] - cyto[0] * (paramsi[6] / paramsi[7])));
  __m128 p02 = _mm_set1_ps(paramsi[2]);
  __m128 p03 = _mm_set1_ps(paramsi[3]);
  __m128 p04 = _mm_set1_ps(paramsi[4]);

  paramsi += nparams;

  __m128 p10 = _mm_set1_ps(paramsi[0]);
  __m128 c10 = _mm_set1_ps(paramsi[1]*(paramsi[5] - cyto[1] * (paramsi[6] / paramsi[7])));
  __m128 p12 = _mm_set1_ps(paramsi[2]);
  __m128 p13 = _mm_set1_ps(paramsi[3]);
  __m128 p14 = _mm_set1_ps(paramsi[4]);

  for (register size_t i = 0; i < npts; i+=4) {
    const size_t j = i + npts;
    
    __m128 x_i = _mm_load_ps(&x_prev[i]);
    __m128 x_j = _mm_load_ps(&x_prev[j]);
    __m128 dx_i = _mm_load_ps(&dx[i]);
    __m128 dx_j = _mm_load_ps(&dx[j]);
    __m128 ddx_i = _mm_load_ps(&ddx[i]);
    __m128 ddx_j = _mm_load_ps(&ddx[j]);
    __m128 temp_i;
    __m128 temp_j;

    temp_i = FASTPOWS(x_j, p04);
    temp_i = _mm_mul_ps(temp_i, p03);
    temp_i = _mm_mul_ps(temp_i, x_i);


    temp_j = FASTPOWS(x_i, p14);
    temp_j = _mm_mul_ps(temp_j, p13);
    temp_j = _mm_mul_ps(temp_j, x_j);

    x_i = _mm_mul_ps(x_i, p02);
    x_j = _mm_mul_ps(x_j, p12);

    ddx_i = _mm_mul_ps(ddx_i, p00);
    ddx_i = _mm_add_ps(ddx_i, c00);


    ddx_i = _mm_sub_ps(ddx_i, x_i);
    ddx_i = _mm_sub_ps(ddx_i, temp_i);
    ddx_i = _mm_sub_ps(ddx_i, dx_i);

    _mm_store_ps(&fx[i], ddx_i);

    ddx_j = _mm_mul_ps(ddx_j, p10);
    ddx_j = _mm_add_ps(ddx_j, c10);
    ddx_j = _mm_sub_ps(ddx_j, x_j);
    ddx_j = _mm_sub_ps(ddx_j, temp_j);
    ddx_j = _mm_sub_ps(ddx_j, dx_j);

    _mm_store_ps(&fx[j], ddx_j);
  }

  return;
}

static void goehring_step(float t, float *x, float *fx) {

  bilinear_flow(flow, current_flow, t * t_flow, npts, nflow);
  finite_difference(x, dx, ddx, current_flow, x_step, npops, npts);
  trapezoidal_integral(x, x_step, cyto, npts, npops);
  reactions(fx, x, dx, ddx, cyto, params, npts, npops, nparams);

/*  if (t>0 && t<20){
    fprintf(fid, "-------------------%.12f---------------\n", t);
    for (int i=0; i<npts; i++) {
      fprintf(fid, "%e: %e %e: %e\n", x[i], dx[i], ddx[i], fx[i]);
      fprintf(fid, "%e: %e %e: %e\n", x[i+npts], dx[i+npts], ddx[i+npts], fx[i+npts]);
    }
  }*/
}

/* x0, dt, tmax, */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  float *time, *x, *fx, *x0, *tmp, *simulation, *curr_sim, *timing;
  float tmax, output_rate, dt, current_dt, t = 0, clf_const, inv_x_step, 
        time_count = 0;
  float relerr, abserr;
  int ntimes, count = 1, page_size = 500, i, j, niter, ntotal, ndata;
  int flag, prev_flag, data_size[2];
  bool saved = false;
  mwSize m, n;

  if (nrhs != 9) {
    mexErrMsgTxt("Simulation requires 9 input arguments !");
  } else {

    x0 = (float*)mxGetData(prhs[0]);
    npts = mxGetM(prhs[0]);
    npops = mxGetN(prhs[0]);
    ntotal = npts*npops;
    ndata = ntotal*sizeof(float);

    params = (float*)mxGetData(prhs[1]);
    nparams = mxGetM(prhs[1]);
    if (mxGetN(prhs[1]) != npops) {
      mexErrMsgTxt("Parameters have to be defined for every population");
    }

    x_step = (float) mxGetScalar(prhs[2]);
    inv_x_step = 1 / x_step;
    tmax = (float) mxGetScalar(prhs[3]);
    dt = (float) mxGetScalar(prhs[4]);
    output_rate = (float) mxGetScalar(prhs[5]);

    if (dt == 0) {
      mexErrMsgTxt("A time step has to be defined");
    }

    flow = (float*)mxGetData(prhs[6]);
    nflow = mxGetN(prhs[6]);
    if (mxGetM(prhs[6]) != npts) {
      mexErrMsgTxt("Flow has to be defined at every position of the lattice !");
    }

    t_flow = (float) (1 / mxGetScalar(prhs[7]));
    ntimes = ceil(tmax / output_rate) + 1;

    niter = (int)(mxGetScalar(prhs[8]));

    if ((x = (float*) _mm_malloc(ndata, 16)) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }

    if ((fx = (float*) _mm_malloc(ndata, 16)) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((dx = (float*) _mm_malloc(ndata, 16)) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((ddx = (float*) _mm_malloc(ndata, 16)) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((current_flow = (float *)mxCalloc(npts, sizeof(float))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((cyto = (float *)mxCalloc(npops, sizeof(float))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((simulation = (float *)mxCalloc(ntotal*ntimes, sizeof(float))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    if ((timing = (float *)mxCalloc(ntimes, sizeof(float))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }

    memcpy(simulation, x0, ndata);
    curr_sim = simulation + ntotal;
    memcpy(x, x0, ndata);
    timing[0] = 0;

    abserr = sqrt(r8_epsilon());
    relerr = sqrt(r8_epsilon());

    prev_flag = -1;

    fid = fopen("sp_run.txt", "w");

    //niter = 200;
    for (j = 0; j < niter; ++j) {
      saved = false;
      current_dt = dt;

      flag = r4_rkf45(goehring_step, ntotal, x, fx, &t, tmax, &relerr, abserr, prev_flag);

      if (flag == 666) {
        char buffer[25];
        sprintf(buffer, "Error due to flag %d\n", prev_flag);
        mexWarnMsgTxt(buffer);
        break;
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

          if ((simulation = (float *)mxRealloc(simulation, ntimes * ndata)) == NULL) {
            mexErrMsgTxt("Memory allocation failed !");
          }
          if ((timing = (float *)mxRealloc(timing, ntimes * sizeof(float))) == NULL) {
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

    fflush(fid);
    fclose(fid);

    data_size[0] = (int)npops*npts;
    data_size[1] = (int)count;
    plhs[0] = mxCreateNumericArray(2, data_size, mxSINGLE_CLASS, mxREAL);
    tmp = (float*)mxGetData(plhs[0]);
    memcpy(tmp, simulation, count*ndata);
    data_size[0] = 1;
    plhs[1] = mxCreateNumericArray(2, data_size, mxSINGLE_CLASS, mxREAL);
    tmp = (float*)mxGetData(plhs[1]);
    memcpy(tmp, timing, count*sizeof(float));

    _mm_free((void*) x);
    _mm_free((void*) dx);
    _mm_free((void*) ddx);
    _mm_free((void*) fx);
    mxFree(cyto);
    mxFree(current_flow);
    mxFree(simulation);
    mxFree(timing);
  }

  return;
}
