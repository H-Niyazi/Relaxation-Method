#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "lattice.h"
#include "charge_cage.h"
#include "plot.h"

// ****************************************************************************
// ************************  charge_cage function  ****************************
// ****************************************************************************
// This method is made based on the Gauss method for the case were we have a
// rectangle of dimensions nx/2 * ny/2 centered at the center of the lattice
// over which the charge e is distributed.
void charge_cage(double * ptr_old, int nx, int ny, int iter, double e)
{

  printf("Charge-cage method was implemented for %d iterations.\n", iter);

  // new pointer to the lattice sites
  double * ptr_new;

  // pointer for saving the data for plotting
  double * ptr_plot;

  // allocating memory for the new pointer and plotting pointer
  ptr_new  = (double *) malloc( nx * ny * sizeof(double) );
  ptr_plot = (double *) malloc( nx * ny * sizeof(double) * ( iter + 1 ) );

  // initializing the first nx * ny elememnts of ptr_plot
  for (int j = 0; j < nx * ny; j++)
    ptr_plot[j] = ptr_old[j];

  // cage perimeter and the charge density over it
  int length = nx + ny - 4;
  double q = e / length;


  for (int i = 1; i < iter + 1; i++)
  {

    for (int j = 0; j < nx * ny; j++)
      ptr_new[j] = ptr_old[j];

    // change only the sites inside the lattice and not the boundaries
    for (int x = 1; x < nx-1; x++)
    {
      for (int y = 1; y < ny-1; y++)
      {
        int site    = idx(x,y,nx);
        int site_xp = idx(x+1,y,nx);
        int site_xm = idx(x-1,y,nx);
        int site_yp = idx(x,y+1,nx);
        int site_ym = idx(x,y-1,nx);

        ptr_new[site] = 0.25 * ( ptr_new[site_xm] + ptr_old[site_xp]
                               + ptr_new[site_ym] + ptr_old[site_yp]);

        // adding the charge-cage term
        if (   (x == nx/4 && y > ny/4 - 1 && y < 3*ny/4)
            || (x == 3*nx/4 - 1 && y > ny/4 - 1 && y < 3*ny/4)
            || (y == ny/4 && x > nx/4 - 1 && x < 3*nx/4)
            || (y == 3*ny/4 - 1 && x > nx/4 - 1 && x < 3*nx/4) )
          ptr_new[site] += 3.14159265359 * q;
      }
    }

    for (int j = 0; j < nx * ny; j++)
      ptr_old[j] = ptr_new[j];

    // adding more elememnts to ptr_plot
    for (int j = 0; j < nx * ny; j++)
      ptr_plot [ j + i * nx * ny ] = ptr_old [ j ];

  }

  // making the plot of the data
  plot_gif (ptr_plot, iter, nx, ny, "Charge-cage Method");

  // deallocating the memory defined inside this function
  free (ptr_new);
  free (ptr_plot);

}

// ****************************************************************************
// *********************  charge_cage_iteration function  *********************
// ****************************************************************************
int charge_cage_iteration(double * ptr, int nx, int ny, double e, double prc, int iter_max)
{

  int iter = 0;

  // sum of the potential on all the inner sites of the lattice
  // initialisation is chosen so that it fails the convergence
  // criterion which is : | sum_new - sum_old | < nx * ny * prc
  double sum_old = 0;
  double sum_new = 2 * nx * ny * prc;

  // old and new pointers to the lattice sites
  double * ptr_old;
  double * ptr_new;

  // allocating memory for the old and new pointers
  ptr_old  = (double *) malloc( nx * ny * sizeof(double) );
  ptr_new  = (double *) malloc( nx * ny * sizeof(double) );

  // initializing the old pointer
  for (int i = 0; i < nx * ny; i++)
    ptr_old[i] = ptr[i];

  // cage perimeter and the charge density over it
  int length = nx + ny - 4;
  double q = e / length;


  // iteration finder loop
  while ( iter < iter_max && ( sum_new - sum_old > nx * ny * prc ||
                               sum_old - sum_new > nx * ny * prc   ) )
  {
    for (int i = 0; i < nx * ny; i++)
      ptr_new[i] = ptr_old[i];

    // for checking the convergence criterion
    sum_old = 0;
    sum_new = 0;

    // change only the sites inside the lattice and not the boundaries
    for (int x = 1; x < nx-1; x++)
    {
      for (int y = 1; y < ny-1; y++)
      {
        int site    = idx(x,y,nx);
        int site_xp = idx(x+1,y,nx);
        int site_xm = idx(x-1,y,nx);
        int site_yp = idx(x,y+1,nx);
        int site_ym = idx(x,y-1,nx);

        ptr_new[site] = 0.25 * ( ptr_new[site_xm] + ptr_old[site_xp]
                               + ptr_new[site_ym] + ptr_old[site_yp]);

        // adding the charge-cage term
        if (   (x == nx/4 && y > ny/4 - 1 && y < 3*ny/4)
            || (x == 3*nx/4 - 1 && y > ny/4 - 1 && y < 3*ny/4)
            || (y == ny/4 && x > nx/4 - 1 && x < 3*nx/4)
            || (y == 3*ny/4 - 1 && x > nx/4 - 1 && x < 3*nx/4) )
          ptr_new[site] += 3.14159265359 * q;

        sum_old = sum_old + ptr_old[site];
        sum_new = sum_new + ptr_new[site];
      }
    }

    for (int i = 0; i < nx * ny; i++)
      ptr_old[i] = ptr_new[i];

    iter = iter + 1;
  }

  if (iter == iter_max)
    printf("Couldn't reach the desired accuracy after %d iterations for the Charge-cage method.\n", iter);

  return iter;

  // deallocating the memory defined inside this function
  free (ptr_old);
  free (ptr_new);

}
