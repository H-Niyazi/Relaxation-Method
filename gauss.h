#ifndef GAUSS_H
#define GAUSS_H

int gauss_iteration(double * ptr, int nx, int ny, int src[2], double e, double prc, int iter_max);
void gauss(double * ptr_old, int nx, int ny, int iter, int src[2], double e);

#endif
