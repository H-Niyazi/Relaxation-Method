#ifndef OVER_RELAXATION_H
#define OVER_RELAXATION_H

int over_relaxation_iteration(double * ptr, int nx, int ny, int src[2], double e, double omega, double prc, int iter_max);
void over_relaxation(double * ptr_old, int nx, int ny, int iter, int src[2], double e, double omega);
double over_relaxation_omega(double * ptr, int nx, int ny, int src[2], double e, double delta_omega, double prc, int iter_max);

#endif
