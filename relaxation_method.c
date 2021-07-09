#include <stdio.h>
#include <stdlib.h>
#include "lattice.h"
#include "jacob.h"
#include "gauss.h"
#include "over_relaxation.h"
#include "charge_cage.h"
#include "plot.h"


int main()
{

  // lattice dimensions
  int nx;
  int ny;
  printf("\nEnter the dimensions of the lattice : \n");
  printf("nx = ");
  scanf("%d", &nx);
  printf("ny = ");
  scanf("%d", &ny);

  // boundary conditions
  double bc_l = 1;  // left  boundary
  double bc_r = 0;  // right boundary
  double bc_u = 1;  // upper boundary
  double bc_d = 0;  // down  boundary

  // delta source position and the scaled magnitude
  int source[2] = { nx / 2, ny / 2 };
  double e;
  printf("Enter the magnitude of the charge : \n");
  scanf("%lf", &e);

  // number of sweeps through the lattice and its maximum
  int iter;
  int iter_max = 10000;

  // precision for the convergence criterion
  double prc = 0.0001;
  printf("\nconvergence precision =  %f \n\n", prc);

  // the best omega and the stepsize for seraching for it for the
  // over-relaxation method
  double delta_omega = 0.001;
  double omega_crit;

  // lattice pointer and its allocated memory
  double * latPtr;
  latPtr = (double *) malloc ( nx * ny * sizeof(double) );


  // finding the number of iterations for the Jacobi method and implementing it
  initializer (nx, ny, latPtr, bc_l, bc_r, bc_u, bc_d);
  iter = jacob_iteration(latPtr, nx, ny, source, e, prc, iter_max);
  jacob(latPtr, nx, ny, iter, source, e);


  // finding the number of iterations for the Gauss method and implementing it
  initializer (nx, ny, latPtr, bc_l, bc_r, bc_u, bc_d);
  iter = gauss_iteration(latPtr, nx, ny, source, e, prc, iter_max);
  gauss(latPtr, nx, ny, iter, source, e);


  // finding the number of iterations and the best omega for the over-relaxation
  // method and implementing it
  initializer (nx, ny, latPtr, bc_l, bc_r, bc_u, bc_d);
  omega_crit = over_relaxation_omega(latPtr, nx, ny, source, e, delta_omega, prc, iter_max);
  iter = over_relaxation_iteration(latPtr, nx, ny, source, e, omega_crit, prc, iter_max);
  over_relaxation(latPtr, nx, ny, iter, source, e, omega_crit);

  // finding the number of iterations for the charge-cage method and
  // implementing it. Here In order to get a more symmetric plot, instead of
  // using the above-defined boundary conditions, I use 0 for all the potentail
  // on all the four boundaries.
  initializer (nx, ny, latPtr, 0, 0, 0, 0);
  iter = charge_cage_iteration(latPtr, nx, ny, e, prc, iter_max);
  charge_cage(latPtr, nx, ny, iter, e);


  // deallocate the memory
  free(latPtr);

}
