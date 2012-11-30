#include <xmmintrin.h> //SSE
#include <emmintrin.h> //SSE2
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "mex.h"
#include "rkf45.h"

#define POS(x) ((x) > 0 ? 1 : 0)
#define NEG(x) ((x) < 0 ? 1 : 0)

#define FASTPOWS(x, y) exp2f4(_mm_mul_ps(log2f4(x), (y)))

#define EXP_POLY_DEGREE 3
#define LOG_POLY_DEGREE 5

#define POLY0(x, c0) _mm_set1_ps(c0)
#define POLY1(x, c0, c1) _mm_add_ps(_mm_mul_ps(POLY0(x, c1), x), _mm_set1_ps(c0))
#define POLY2(x, c0, c1, c2) _mm_add_ps(_mm_mul_ps(POLY1(x, c1, c2), x), _mm_set1_ps(c0))
#define POLY3(x, c0, c1, c2, c3) _mm_add_ps(_mm_mul_ps(POLY2(x, c1, c2, c3), x), _mm_set1_ps(c0))
#define POLY4(x, c0, c1, c2, c3, c4) _mm_add_ps(_mm_mul_ps(POLY3(x, c1, c2, c3, c4), x), _mm_set1_ps(c0))
#define POLY5(x, c0, c1, c2, c3, c4, c5) _mm_add_ps(_mm_mul_ps(POLY4(x, c1, c2, c3, c4, c5), x), _mm_set1_ps(c0))

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
float flow_step, x_step;
int npts, npops, nflow, nparams, ncalls;

static void finite_difference(float *t_x, float *t_dx, float *t_ddx, float *t_flow, float t_h, int t_npop, int t_npts) {

  float dnorm, ddnorm;
  register size_t i;

  dnorm = 1 / t_h;
  ddnorm = 1 / (12*t_h*t_h);

  __m128 x_i1, mid1, tmp1;

  __m128 dn = _mm_set1_ps(dnorm);
  __m128 ddn = _mm_set1_ps(ddnorm);
  __m128 z = _mm_set1_ps(0.0f);
  __m128 c1 = _mm_set1_ps(16.0f);
  __m128 c0 = _mm_set1_ps(30.0f);
  __m128 u1 = _mm_set1_ps(4.0f);
  __m128 u0 = _mm_set1_ps(3.0f);

  __m128 x_i = _mm_load_ps(&t_x[0]);
  __m128 x_i_1 = _mm_shuffle_ps(x_i, x_i, _MM_SHUFFLE(1,0,3,2));
  __m128 mid_1 = _mm_shuffle_ps(x_i_1, x_i, _MM_SHUFFLE(1,0,3,2));

  __m128 f_i = _mm_load_ps(&t_flow[0]);
  __m128 f_i_1 = _mm_shuffle_ps(f_i, f_i, _MM_SHUFFLE(1,0,3,2));
  __m128 midf_1 = _mm_shuffle_ps(f_i_1, f_i, _MM_SHUFFLE(1,0,3,2));

  midf_1 = _mm_mul_ps(mid_1, midf_1);
  f_i_1 = _mm_mul_ps(x_i_1, f_i_1);
  f_i = _mm_mul_ps(x_i, f_i);

  for (i=0; i<((t_npop*t_npts)-4); i+=4) {

    x_i1 = _mm_load_ps(&t_x[i+4]);
    mid1 = _mm_shuffle_ps(x_i, x_i1, _MM_SHUFFLE(1,0,3,2));

    tmp1 = _mm_add_ps(mid_1, mid1);
    tmp1 = _mm_mul_ps(tmp1, c1);
    x_i_1 = _mm_add_ps(x_i_1, x_i1);
    x_i_1 = _mm_sub_ps(tmp1, x_i_1);
    tmp1 = _mm_mul_ps(x_i, c0);
    x_i_1 = _mm_sub_ps(x_i_1, tmp1);
    x_i_1 = _mm_mul_ps(x_i_1, ddn);

    _mm_store_ps(&t_ddx[i], x_i_1);
    
    __m128 f_i1 = _mm_load_ps(&t_flow[i+4]);
    f_i1 = _mm_mul_ps(x_i1, f_i1);
    __m128 midf1 = _mm_shuffle_ps(f_i, f_i1, _MM_SHUFFLE(1,0,3,2));

    __m128 neg = _mm_cmplt_ps(f_i, z);
    __m128 pos = _mm_cmpgt_ps(f_i, z);

    midf_1 = _mm_mul_ps(midf_1, u1);
    tmp1 = _mm_mul_ps(f_i, u0);
    midf_1 = _mm_sub_ps(tmp1, midf_1);
    midf_1 = _mm_add_ps(midf_1, f_i_1);
    midf_1 = _mm_and_ps(midf_1, pos);

    f_i_1 = _mm_add_ps(tmp1, f_i1);
    tmp1 = _mm_mul_ps(midf1, u1);
    f_i_1 = _mm_sub_ps(tmp1, f_i_1);
    f_i_1 = _mm_and_ps(f_i_1, neg);

    tmp1 = _mm_or_ps(midf_1, f_i_1);
    tmp1 = _mm_mul_ps(tmp1, dn);

    _mm_store_ps(&t_dx[i], tmp1);

    midf_1 = midf1;
    f_i_1 = f_i;
    f_i = f_i1;

    x_i_1 = x_i;
    x_i = x_i1;
    mid_1 = mid1;
  }

  x_i1 = _mm_shuffle_ps(x_i, x_i, _MM_SHUFFLE(1,0,3,2));
  mid1 = _mm_shuffle_ps(x_i, x_i1, _MM_SHUFFLE(1,0,3,2));

  tmp1 = _mm_add_ps(mid_1, mid1);
  tmp1 = _mm_mul_ps(tmp1, c1);
  x_i_1 = _mm_add_ps(x_i_1, x_i1);
  x_i_1 = _mm_sub_ps(tmp1, x_i_1);
  tmp1 = _mm_mul_ps(x_i, c0);
  x_i_1 = _mm_sub_ps(x_i_1, tmp1);

  x_i_1 = _mm_mul_ps(x_i_1, ddn);

  _mm_store_ps(&t_ddx[i], x_i_1);

  __m128 f_i1 = _mm_shuffle_ps(f_i, f_i, _MM_SHUFFLE(1,0,3,2));
  __m128 midf1 = _mm_shuffle_ps(f_i, f_i1, _MM_SHUFFLE(1,0,3,2));

  __m128 neg = _mm_cmplt_ps(f_i, z);
  __m128 pos = _mm_cmpgt_ps(f_i, z);

  midf_1 = _mm_mul_ps(midf_1, u1);
  tmp1 = _mm_mul_ps(f_i, u0);
  midf_1 = _mm_sub_ps(tmp1, midf_1);
  midf_1 = _mm_add_ps(midf_1, f_i_1);
  midf_1 = _mm_and_ps(midf_1, pos);

  f_i_1 = _mm_add_ps(tmp1, f_i1);
  tmp1 = _mm_mul_ps(midf1, u1);
  f_i_1 = _mm_sub_ps(tmp1, f_i_1);
  f_i_1 = _mm_and_ps(f_i_1, neg);
  tmp1 = _mm_or_ps(midf_1, f_i_1);
  tmp1 = _mm_mul_ps(tmp1, dn);

  _mm_store_ps(&t_dx[i], tmp1);

  return;
}

static void bilinear_flow(float *t_flow, float *t_current, float t_index, int t_npts, int t_ntimes) {

  float yf, yc;
  int indf, indc;

  if (t_index >= t_ntimes-1) {
    memset(t_current, 0.0f, 2*t_npts*sizeof(float));

    return;
  } else {
    indf = (int)floor(t_index);
    yf = 1 - (t_index - (float)indf);

    if (yf == 0) {
      indc = indf;
      yc = 1;
    } else {
      indc = indf + 1;
      yc = 1 - yf;
    }
  }

  indf = indf*t_npts;
  indc = indc*t_npts;

  __m128 y_f = _mm_set1_ps(yf);
  __m128 y_c = _mm_set1_ps(yc);

  for (register size_t i=0, j=0; i<t_npts; i+=4, j+=8) {
    __m128 f_f = _mm_load_ps(&t_flow[i+indf]);
    __m128 f_c = _mm_load_ps(&t_flow[i+indc]);
    f_f = _mm_mul_ps(f_f, y_f);
    f_c = _mm_mul_ps(f_c, y_c);
    f_f = _mm_add_ps(f_f, f_c);
    f_c = _mm_shuffle_ps(f_f, f_f, _MM_SHUFFLE(3,3,2,2));
    f_f = _mm_shuffle_ps(f_f, f_f, _MM_SHUFFLE(1,1,0,0));

    _mm_store_ps(&t_current[j], f_f);
    _mm_store_ps(&t_current[j+4], f_c);
  }

  return;
}

static void trapezoidal_integral(float *t_x, float t_h, float *t_ints, int t_npts, int t_npops) {

  float tmp_0, tmp_1;
  int i, p;

  tmp_0 = t_x[0];
  tmp_1 = t_x[1];
  for (i=2; i<t_npops*(t_npts-1); i+=2) {
    tmp_0 += 2*t_x[i];
    tmp_1 += 2*t_x[i+1];
  }
  tmp_0 += t_x[i];
  tmp_1 += t_x[i+1];

  t_ints[0] = tmp_0 * (t_h/2);
  t_ints[1] = tmp_1 * (t_h/2);

  return;
}

static void reactions(float *t_fx, float *t_x_prev, float *t_dx, float *t_ddx, float *t_cyto, float *t_params, int t_npts, int t_npops, int t_nparams) {

  float *paramsi;

  paramsi = t_params;

  __m128 p00 = _mm_set1_ps(paramsi[0]);
  __m128 c00 = _mm_set1_ps(paramsi[1]*(paramsi[5] - t_cyto[0] * (paramsi[6] / paramsi[7])));
  __m128 p02 = _mm_set1_ps(paramsi[2]);
  __m128 p03 = _mm_set1_ps(paramsi[3] * paramsi[4 + t_nparams]);
  __m128 p04 = _mm_set1_ps(paramsi[4]);

  paramsi += t_nparams;

  __m128 p10 = _mm_set1_ps(paramsi[0]);
  __m128 c10 = _mm_set1_ps(paramsi[1]*(paramsi[5] - t_cyto[1] * (paramsi[6] / paramsi[7])));
  __m128 p12 = _mm_set1_ps(paramsi[2]);
  __m128 p13 = _mm_set1_ps(paramsi[3] * params[4]);
  __m128 p14 = _mm_set1_ps(paramsi[4]);

  __m128 p0 = _mm_unpacklo_ps(p00, p10);
  __m128 c0 = _mm_unpacklo_ps(c00, c10);
  __m128 p2 = _mm_unpacklo_ps(p02, p12);
  __m128 p3 = _mm_unpacklo_ps(p03, p13);
  __m128 p4 = _mm_unpacklo_ps(p04, p14);

  for (register size_t i = 0; i < t_npops*t_npts; i+=4) {
    
    __m128 x_i = _mm_load_ps(&t_x_prev[i]);
    __m128 dx_i = _mm_load_ps(&t_dx[i]);
    __m128 ddx_i = _mm_load_ps(&t_ddx[i]);
    __m128 temp_i = _mm_shuffle_ps(x_i, x_i, _MM_SHUFFLE(2,3,0,1));

    //----------------------------------------------------------
    // _mm_shufle(a,b,_MM_SHUFFLE(i,j,k,l)) produces a[l] a[k] b[j] [b[i]
    //----------------------------------------------------------

    temp_i = FASTPOWS(temp_i, p4);
    temp_i = _mm_mul_ps(temp_i, p3);
    temp_i = _mm_mul_ps(temp_i, x_i);

    x_i = _mm_mul_ps(x_i, p2);

    ddx_i = _mm_mul_ps(ddx_i, p0);
    ddx_i = _mm_add_ps(ddx_i, c0);

    ddx_i = _mm_sub_ps(ddx_i, x_i);
    ddx_i = _mm_sub_ps(ddx_i, temp_i);
    ddx_i = _mm_sub_ps(ddx_i, dx_i);

    _mm_store_ps(&t_fx[i], ddx_i);
  }

  return;
}

static void goehring_step(float t, float *x, float *fx) {

  ncalls++;

  bilinear_flow(flow, current_flow, t * flow_step, npts, nflow);
  finite_difference(x, dx, ddx, current_flow, x_step, npops, npts);
  trapezoidal_integral(x, x_step, cyto, npts, npops);
  reactions(fx, x, dx, ddx, cyto, params, npts, npops, nparams);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  float *time, *x, *fx, *x0, *tmp, *simulation, *curr_sim, *timing, *tmp_simul;
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

    nflow = mxGetN(prhs[6]);
    if (mxGetM(prhs[6]) != npts) {
      mexErrMsgTxt("Flow has to be defined at every position of the lattice !");
    }

    if ((flow = (float*) _mm_malloc(nflow*npts * sizeof(float), 16)) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }
    memcpy(flow, (float*)mxGetData(prhs[6]), nflow*npts*sizeof(float));

    flow_step = (float) (1 / mxGetScalar(prhs[7]));
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
    if ((current_flow = (float *)mxCalloc(ntotal, sizeof(float))) == NULL) {
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
    
    for (j = 0; j < npts; ++j) {
      x[2*j] = x0[j];
      x[2*j+1] = x0[j+npts];
    }
    //memcpy(x, x0, ndata);

    timing[0] = 0;
    ncalls = 0;

    abserr = sqrt(r8_epsilon());
    relerr = sqrt(r8_epsilon());

    prev_flag = -1;
    for (j = 0; j < niter; ++j) {
      saved = false;
      current_dt = dt;

      flag = r4_rkf45(goehring_step, ntotal, x, fx, &t, tmax, &relerr, abserr, prev_flag);
      //mexPrintf("%d / %d\n", ncalls, niter);

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

      if (t >= tmax || ncalls >= niter) {
        break;
      }
    }
    if (!saved) {
      memcpy(curr_sim, x, ndata);
      timing[count] = t;
      count++;
    }
    if (j >= niter || ncalls >= niter) {
      mexWarnMsgTxt("Simulation reached the iteration limit !");
    }

    data_size[0] = (int)npops*npts;
    data_size[1] = (int)count;
    plhs[0] = mxCreateNumericArray(2, data_size, mxSINGLE_CLASS, mxREAL);
    tmp = (float*)mxGetData(plhs[0]);

    memcpy(tmp, simulation, ndata);

    tmp_simul = simulation;
    for (i=1;i<count;i++){
      tmp += ntotal;
      tmp_simul += ntotal;
      for (j=0; j<npts; j++) {
        tmp[j] = tmp_simul[2*j];
        tmp[j+npts] = tmp_simul[2*j+1];
      }
    }

    //memcpy(tmp, simulation, count*ndata);
    data_size[0] = 1;
    plhs[1] = mxCreateNumericArray(2, data_size, mxSINGLE_CLASS, mxREAL);
    tmp = (float*)mxGetData(plhs[1]);
    memcpy(tmp, timing, count*sizeof(float));

    _mm_free((void*) flow);
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
