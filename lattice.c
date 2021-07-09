#include "lattice.h"


// ****************************************************************************
// ***************************  idx function  *********************************
// ****************************************************************************
// idx gives the index of the lattice sites and starts from zero and goes up to
// nx * ny - 1 . It starts from the left down corner (0,0) and goes right to the
// (0,ny-1), then from (0,1) and goes right to the (1,ny-1) and so on. The order
// in which we move can be schematically shown as this:
// ... -------->
// 3   -------->
// 2   -------->
// 1   -------->
// x takes the values from 0 to nx-1 and y takes th evalues from 0 to ny-1
int idx ( int x, int y, int nx)
{
  return x + y * nx;
}

// ****************************************************************************
// ***************************  posx function  ********************************
// ****************************************************************************
// posx gives the x component of each of the lattice sites and gives a number
// from 0 to nx-1 . We have to provide it with the lattice site index which is
// a number from 0 to nx*ny-1 and is obtained by the way defined in the comment
// for the idx function
int posx ( int index, int nx )
{
  return index % nx;
}

// ****************************************************************************
// ***************************  posy function  ********************************
// ****************************************************************************
// posy gives the y component of each of the lattice sites and gives a number
// from 0 to ny-1 . We have to provide it with the lattice site index which is
// a number from 0 to nx*ny-1 and is obtained by the way defined in the comment
// for the idx function
int posy ( int index, int nx )
{
  return index / nx;
}

// ****************************************************************************
// *************************  initializer function  ***************************
// ****************************************************************************
// This function initializes the boundary sites of the lattice to the boundary
// condition values and other sites to zero
void initializer(int nx, int ny, double * latPtr, double bc_l, double bc_r, double bc_u, double bc_d)
{
  // all the sites
  for (int site = 0; site < nx*ny; site++)
    latPtr[site] = 0;

  // boundary sites
  for (int site = 0; site < nx*ny; site++)
  {
    if( posy(site, nx) == 0 )     latPtr[site] = bc_d;
    if( posx(site, nx) == nx-1 )  latPtr[site] = bc_r;
    if( posy(site, nx) == ny-1 )  latPtr[site] = bc_u;
    if( posx(site, nx) == 0 )     latPtr[site] = bc_l;
  }

}
