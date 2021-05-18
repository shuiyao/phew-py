#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include "loadhdf5.h"
#include "gadgetdefs.h"

/* struct particle_data *P; */
/* struct gadget_dump h; */

#define NUMOFNODES_X 256
#define NUMOFNODES_Y 256
#define HUBBLEPARAM 0.7
#define OMEGABARYON 0.045
#define XMIN -1.5
#define XMAX 7.0
#define YMIN 3.0
#define YMAX 8.0

typedef struct NodesStruct {
  float x;
  float y;
  int z;
  double mass;
  double mwind;
} NodesStruct;
NodesStruct *nodes_array;

int ncells_x;
int ncells_y;
char infilename[200], outfilename[200];

int LoadGrid(int nx, int ny)
{
  int i, ix, iy;
  float dx, dy, x;
  int ncells;
  ncells = nx * ny;
  nodes_array = (struct NodesStruct *)
    malloc(ncells*sizeof(*nodes_array));
  dx = (XMAX - XMIN) / nx;
  dy = (YMAX - YMIN) / ny;
  for(ix=0;ix<nx;ix++) {
    x = XMIN + dx * ix;
    for (iy=0;iy<ny;iy++){
      i = ix * ny + iy;
      nodes_array[i].x = x;
      nodes_array[i].y = YMIN + dy * iy;
      nodes_array[i].z = 0;
      nodes_array[i].mass = 0.0;
      nodes_array[i].mwind = 0.0;      
    }
  }
  printf("Grid Loaded.\n");
  return 1;
}

void get_filenames(char *snapbase, int snapnum)
{
  char snapstr[4];
  get_snap_string(snapnum, snapstr);
  sprintf(outfilename, "mrhot_phew_%s", snapstr);
  sprintf(infilename, "%s_%s", snapbase, snapstr);
}

int GridCount(float x, float y, double mass, double wmass)
{
  int ix, iy, i;
  ix = (int)((x - XMIN) / (XMAX - XMIN) * ncells_x);
  if (ix < ncells_x && ix > 0){
    iy = (int)((y - YMIN) / (YMAX - YMIN) * ncells_y);
    if (iy < ncells_y && iy > 0){
      i = ix * ncells_y + iy;
      nodes_array[i].z ++;
      nodes_array[i].mass += mass;
      nodes_array[i].mwind += wmass;      
      return 0;
    }
  }
  return 1;
}

void ReadData(char *filename)
{
  int i, Nout = 0, Ntot = 0;
  float LogRho, LogT;
  double Mass, WindMass;
  double MeanWeight;

  load_hdf5(filename);
  cosmounits();
  
  for(i=0;i<h.npart[0];i++) 
    {
      LogRho = P[i].Rho * UNIT_M / pow(UNIT_L, 3) / unit_Density;
      LogRho = log10(LogRho * HUBBLEPARAM * HUBBLEPARAM / OMEGABARYON);
      MeanWeight = (1 + 4 * XHE) / (1 + P[i].Ne + XHE);
      LogT = P[i].Temp * unit_Temp;
      LogT *= GAMMA_MINUS1 * PROTONMASS / BOLTZMANN * MeanWeight;
      LogT = log10(LogT);
      Mass = P[i].Mass;
#ifdef PHEW_EXTRA_OUTPUT      
      if(P[i].Mcloud > 0) WindMass = P[i].Mass;
      else WindMass = P[i].WindMass;
#endif      
      Nout += GridCount(LogRho, LogT, Mass, WindMass);
      Ntot ++;
    }
  printf("Data Reading Done, Nout/Ntot = %d/%d\n", Nout, Ntot);
}

void WriteGrid(char *filename)
{
  FILE *outfile;
  int ix, iy, i;
  outfile = fopen(filename, "w");
  fprintf(outfile, "%d %d\n", ncells_x, ncells_y);
  for(ix=0;ix<ncells_x;ix++) {fprintf(outfile, "%g\n", nodes_array[ix*ncells_y].x);}
  for(iy=0;iy<ncells_y;iy++) {fprintf(outfile, "%g\n", nodes_array[iy].y);}
  for(ix=0;ix<ncells_x;ix++) {
    for (iy=0;iy<ncells_y;iy++){
      i = ix * ncells_y + iy;
      fprintf(outfile, "%d %g %g\n", nodes_array[i].z,
	      nodes_array[i].mass,
	      (nodes_array[i].mass) ? nodes_array[i].mwind / nodes_array[i].mass : 0.0);
    }
  }
  fclose(outfile);
}

int main(int argc, char **argv)
{
  /* Calling Sequence: */
  /* rhot_hist snapbase 00 33 256 256 */
  char snapbase[200];
  int snapnum;  
  strcpy(snapbase, argv[1]);
  snapnum = atoi(argv[2]);
  if (argc == 5)
    { ncells_x = atoi(argv[3]);
      ncells_y = atoi(argv[4]);
    }
  else
    { ncells_x = NUMOFNODES_X;
      ncells_y = NUMOFNODES_Y;
    }
  get_filenames(snapbase, snapnum);
  LoadGrid(ncells_x, ncells_y);
  ReadData(infilename);
  WriteGrid(outfilename);
  free(nodes_array);
  return 1;
}
