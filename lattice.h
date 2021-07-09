#ifndef LATTICE_H
#define LATTICE_H

int idx ( int x, int y, int nx);
int posx (int index, int nx);
int posy (int index, int nx);

void initializer(int nx, int ny, double * latPtr, double bc_l, double bc_r, double bc_u, double bc_d);

#endif
