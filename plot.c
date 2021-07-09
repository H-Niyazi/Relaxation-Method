#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>     // required for chdir()
#include <sys/types.h>  // required for stat()
#include <sys/stat.h>   // required for stat()
#include "lattice.h"
#include "plot.h"

// this code is written with the help from these webpages:
// http://gnuplot-surprising.blogspot.com/2011/09/creating-gif-animation-using-gnuplot.html

// ****************************************************************************
// ***************************  plot_gif function  ****************************
// ****************************************************************************
void plot_gif ( double * ptr_plot, int iter, int nx, int ny, char * title)
{

  // making the data directory
  struct stat st = {0};
  const char * dir = title;
  if (stat(dir, &st) == -1)   // check if the data directory already exists
    mkdir(dir, 0700);


  // making the data files inside the data directory
  for (int j = 0; j < iter + 1; j++)
  {
    char name[100];
    sprintf (name, "%s/data%d.txt", dir, j);
    FILE * data = fopen ( name, "w" );

    for (int k = 0; k < nx; k++)
    {
      for (int l = 0; l < ny; l++)
      {
        fprintf (data, "%d\t%d\t%f\n", k, l, ptr_plot[k + l*nx + j*nx*ny]);
      }
      // this newline is necessary for pm3d in gnuplot to work
      fprintf (data, "\n");
    }

    fclose(data);
  }


  // making the file to be loaded to gnuplot inside the data directory
  char animate_name[50];
  sprintf(animate_name, "%s/animate.plt", dir);
  FILE * animate = fopen ( animate_name, "w" );

  fprintf(animate, "splot sprintf(\"data%%i.txt\",i) using 1:2:3 with pm3d title sprintf(\"Iteration = %%i\",i)\n");
  fprintf(animate, "i = i + 1\n");
  fprintf(animate, "if (i < n) reread"); // n is declared iter + 1 later


  // changing the working directory to the data directory
  int chdir(const char * path);
  if (chdir(dir) == 0)
    printf("New directory made successfully for saving the plots.\n");
  else
    printf("Could not change directory.\n");


  // defining the pipe to gnuplot
  FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");

  // gnuplot commands
  fprintf(gnuplotPipe, "reset \n");
  fprintf(gnuplotPipe, "set xlabel \"x\" \n");
  fprintf(gnuplotPipe, "set ylabel \"y\" \n");
  fprintf(gnuplotPipe, "set zlabel \"V\" \n");

  fprintf(gnuplotPipe, "set xr [0:%d] \n", nx);
  fprintf(gnuplotPipe, "set yr [0:%d] \n", ny);
  fprintf(gnuplotPipe, "set title \"%s\" \n", title);
  fprintf(gnuplotPipe, "set contour\n");

  fprintf(gnuplotPipe, "set term gif animate \n");
  fprintf(gnuplotPipe, "set output \"0 %s.gif\" \n", title);
  fprintf(gnuplotPipe, "n=%d \n", iter + 1);  // total number of images in gif
  fprintf(gnuplotPipe, "i=0 \n");                 // iteration over the images
  fprintf(gnuplotPipe, "load \"animate.plt\" \n");     // load plotting script
  fprintf(gnuplotPipe, "set output \n");

  fprintf (gnuplotPipe, "system('rm *.txt') \n");   // remove plotting script
  fprintf (gnuplotPipe, "system('rm *.plt') \n");   // remove data files

  fprintf (gnuplotPipe, "exit \n");   // exit gnuplot


  // changing the working directory to the main directory
  int chdir(const char * back);
  char back[10] = "..";
  if (chdir(back) == 0)
    printf("Successfully changed back to the main directory.\n\n");
  else
    printf("Could not change directory.\n\n");

}
