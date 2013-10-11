/*
 *******************  compute quantile of adaptive Neyman test**********
 * based on the article 
 * Fan, J. and Lin, S.K. (1998), Test of Significance when data are curves., 
 * Journal of American Statistical Association, 93, 1007-1021.
 *
 * adapted from the C code from Jianqing Fan:
 * http://www.orfe.princeton.edu/~jqfan/papers/pub/neyman.c
 */
#include <math.h>
#include <stdlib.h>
#include "mex.h"

#define NSIM 100000
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ******************** UNIFORM RANDOM NUMBER GENERATOR: RAN1.c *******
*/
#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

static double ran1(int *idum) {

                    /* returns a uniform random deviate between 0.0 and
1.0.
                       Set idum to any negative value to initialize or
                       reinitialize the sequence  */

    static long ix1,ix2,ix3;
    static double r[98];
    double temp;
    static int iff=0;
    int j;

    if (*idum < 0 || iff == 0) {
      iff=1;
      ix1=(IC1-(*idum)) % M1;
      ix1=(IA1*ix1+IC1) % M1;
      ix2=ix1 % M2;
      ix1=(IA1*ix1+IC1) % M1;
      ix3=ix1 % M3;
      for (j=1;j<=97;j++) {
        ix1=(IA1*ix1+IC1) % M1;
        ix2=(IA2*ix2+IC2) % M2;
        r[j]=(ix1+ix2*RM2)*RM1;
      }
      *idum=1;
    }
    ix1=(IA1*ix1+IC1) % M1;
    ix2=(IA2*ix2+IC2) % M2;
    ix3=(IA3*ix3+IC3) % M3;

    j = 1 + ((97*ix3)/M3);

    if (j > 97 || j < 1) {
      mexErrMsgTxt("Error in the random number generator !");
    }

    temp=r[j];
    r[j]=(ix1+ix2*RM2)*RM1;

    return temp;
}


/* ************************ADPATIVE NEYMAN.C*********************************
Compute Adaptive Neyman test statisti
x: n by 1 vector, is double.  
************************************************************************** */
static double aneyman(double *x, int n) {

  int i, mopt;
  double sum, teststat, valueopt, valuenew;
  double logn, constant;

  constant = sqrt(2/M_PI);          /* sqrt(2/pi) */

  /*  Computing Adaptive Neyman statistic  */
  mopt = 1;
  sum = x[0]*x[0];
  valueopt = sum - 1.0;

  for(i=1; i < n; i++) {
    sum += x[i]*x[i];
    valuenew = (sum-(i+1))/sqrt((double) (i+1));

    if(valuenew > valueopt) {
      mopt = i+1;
      valueopt = valuenew;
    }
  }

  teststat = valueopt/sqrt(2.0);
  logn = log( log ( ( double) n ) );
  teststat = sqrt(2*logn)*teststat - (2*logn + 0.5*log(logn)-1.26);

  return teststat;
}


/* ***************** NORMAL RANDOM GENERATOR: GASDEV.C *********** */
static double gasdev(int* idum) {
                   /* returns a normally distributed deviate with mean 0
                     and variance 1, using ran1(idum) as the source of uniform
                      deviates    */
  static int iset=0;
  static double gset;
  double fac,r,v1,v2;
  double ran1();

  if (iset == 0)  {
    do {
      v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      r=v1*v1+v2*v2;
    } while (r >= 1.0 || r == 0.0);

    fac= sqrt(-2.0*log(r)/r);
    gset=v1*fac;
    iset=1;

    return v2*fac;
  }  else {
    iset=0;
    return gset;
  }
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  int n, i, nsim, idum;
  double teststat, neyman[NSIM], *x, sum = 0.0, *pval;

  if (nrhs != 2) {
    mexErrMsgTxt("This function accepts only 2 input arguments (mhat and Tan) !");
  } else {

    n = (int)(mxGetScalar(prhs[0]));
    teststat = mxGetScalar(prhs[1]);

    if ((x = mxCalloc(n, sizeof(double))) == NULL) {
      mexErrMsgTxt("Memory allocation failed !");
    }

    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    pval = mxGetPr(plhs[0]);

    idum = -rand();

    for(nsim=0; nsim < NSIM; nsim++) {      /* simulations */
      for(i=0; i < n; i++) {
        x[i] = gasdev(&idum);
      }
      neyman[nsim]  = aneyman(x, n);
    }

    for(i=0; i < NSIM; i++) {
      if(neyman[i] > teststat) {
        sum += 1.0;
      }
    }
    pval[0] = sum/NSIM;           /*  P-values */
  }

  return;
}
