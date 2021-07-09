#ifndef CHARGE_CAGE_H
#define CHARGE_CAGE_H

int charge_cage_iteration(double * ptr, int nx, int ny, double e, double prc, int iter_max);
void charge_cage(double * ptr_old, int nx, int ny, int iter, double e);

#endif
