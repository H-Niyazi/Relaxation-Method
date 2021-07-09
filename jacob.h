#ifndef JACOB_H
#define JACOB_H

int jacob_iteration(double * ptr, int nx, int ny, int src[2], double e, double prc, int iter_max);
void jacob(double * ptr_old, int nx, int ny, int iter, int src[2], double e);

#endif
