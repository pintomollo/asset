#ifndef GAUSS_H
#define GAUSS_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
extern "C" {
#endif

void gaussian_smooth(double *image, int rows, int cols, double sigma);
void make_gaussian_kernel(double sigma, double **kernel, int *windowsize);

#ifdef __cplusplus
}
#endif

#endif
