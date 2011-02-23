#ifndef GAUSS_H
#define GAUSS_H

#ifdef __cplusplus
extern "C" {
#endif

void gaussian_smooth(double *image, int rows, int cols, double sigma);
void make_gaussian_kernel(double sigma, double **kernel, int *windowsize);

#ifdef __cplusplus
}
#endif

#endif
