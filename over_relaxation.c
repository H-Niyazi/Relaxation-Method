#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "lattice.h"
#include "over_relaxation.h"
#include "plot.h"

// ****************************************************************************
// ***************************  over_relaxation function  *********************
// ****************************************************************************
void over_relaxation(double * ptr_old, int nx, int ny, int iter, int src[2], double e, double omega)
{

  printf("Over-relaxation method was implemented for %d iterations and omega = %f.\n", iter, omega);

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

        ptr_new[site] = ( 1 - omega ) * ptr_old[site]
                        + 0.25 * omega * (ptr_new[site_xm] + ptr_old[site_xp]
                                         +ptr_new[site_ym] + ptr_old[site_yp]);

        // adding the delta charge term
        if ( x == src[0] && y == src[1] )
          ptr_new[site] += 3.14159265359 * e * omega;
      }
    }

    for (int j = 0; j < nx * ny; j++)
      ptr_old[j] = ptr_new[j];

    // adding more elememnts to ptr_plot
    for (int j = 0; j < nx * ny; j++)
      ptr_plot [ j + i * nx * ny ] = ptr_old [ j ];

  }

  // making the plot of the data
  plot_gif (ptr_plot, iter, nx, ny, "Over-relaxation Method");

  // deallocating the memory defined inside this function
  free (ptr_new);
  free (ptr_plot);

}

// ****************************************************************************
// ********************  over_relaxation_iteration function  ******************
// ****************************************************************************
int over_relaxation_iteration(double * ptr, int nx, int ny, int src[2], double e, double omega, double prc, int iter_max)
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

        ptr_new[site] = ( 1 - omega ) * ptr_old[site]
                        + 0.25 * omega * (ptr_new[site_xm] + ptr_old[site_xp]
                                         +ptr_new[site_ym] + ptr_old[site_yp]);

        // adding the delta charge term
        if ( x == src[0] && y == src[1] )
          ptr_new[site] += 3.14159265359 * e * omega;

        sum_old = sum_old + ptr_old[site];
        sum_new = sum_new + ptr_new[site];
      }
    }

    for (int i = 0; i < nx * ny; i++)
      ptr_old[i] = ptr_new[i];

    iter = iter + 1;
  }

  if (iter == iter_max)
    printf("Couldn't reach the desired accuracy after %d iterations for the over-relaxation method.\n", iter);

  return iter;

  // deallocating the memory defined inside this function
  free (ptr_old);
  free (ptr_new);

}

// ****************************************************************************
// ************************  over_relaxation_omega function  ******************
// ****************************************************************************
// This function finds the omega for which the over-relaxation method has the
// fastest convergence.
double over_relaxation_omega(double * ptr, int nx, int ny, int src[2], double e, double delta_omega, double prc, int iter_max)
{
  double omega = 1.0;
  double omega_crit;
  int iter_min = iter_max;

  while ( omega < 2.0 + delta_omega )
  {
    int iter = over_relaxation_iteration(ptr, nx, ny, src, e, omega, prc, iter_max);
    if (iter < iter_min)
    {
      iter_min = iter;
      omega_crit = omega;
    }
    omega = omega + delta_omega;
  }

  return omega_crit;
}
