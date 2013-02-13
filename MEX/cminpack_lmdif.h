#include <math.h>
#include <float.h>

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503131e-16
#endif
#ifndef DBL_MIN
#define DBL_MIN 2.2250738585072014e-308
#endif
#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623157e+308
#endif

#define double_dwarf (1.4916681462400413e-154*0.9)
#define double_giant (1.3407807929942596e+154*0.9)

#define p0001 1e-4
#define p001 .001
#define p05 .05
#define p1 .1
#define p25 .25
#define p5 .5
#define p75 .75

int cminpack_lmdif(int (*func)(void *d, int n, int m, double *p, double *hx, int iflag), 
    void *d, int m, int n, double *x, double *fvec, double tol, int *iwa, 
    double *wa, int lwa);

int lmdif(int (*func)(void *d, int n, int m, double *p, double *hx, int iflag), 
    void *d, int m, int n, double *x, 
    double *fvec, double ftol, double xtol, double
    gtol, int maxfev, double epsfcn, double *diag, int
    mode, double factor, int nprint, int *
    nfev, double *fjac, int ldfjac, int *ipvt, double *
    qtf, double *wa1, double *wa2, double *wa3, double *
    wa4);

double enorm(int n, const double *x);

int fdjac2(int (*func)(void *d, int n, int m, double *p, double *hx, int iflag), 
    void *d, int m, int n, double *x, 
    const double *fvec, double *fjac, int ldfjac,
    double epsfcn, double *wa);

void qrfac(int m, int n, double *a, int
    lda, int pivot, int *ipvt, int lipvt, double *rdiag,
    double *acnorm, double *wa);

void lmpar(int n, double *r, int ldr, 
    const int *ipvt, const double *diag, const double *qtb, double delta, 
    double *par, double *x, double *sdiag, double *wa1, 
    double *wa2);

void qrsolv(int n, double *r, int ldr, 
    const int *ipvt, const double *diag, const double *qtb, double *x, 
    double *sdiag, double *wa);
