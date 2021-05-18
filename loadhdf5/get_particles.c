/* 180630 [CRITICAL] Bug Fix. Fix the mass and other properties for stars. */

/* Need Mass Msub, Rsub, WindMass, dR, Mcloud */

#include <stdio.h>
#include "gadgetdefs.h"
#include "loadhdf5.h"
#include <hdf5.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int flag_write_phew = 0;
int flag_write_sph = 0;
int flag_write_halo = 0;
int haloid;
char base_name[MAX_LEN_FILENAME];

#define NSPH 50000
#define OMEGABARYON 0.045

char binname[MAX_LEN_FILENAME];
char grpname[MAX_LEN_FILENAME], gtpname[MAX_LEN_FILENAME];
char sogrpname[MAX_LEN_FILENAME], sogtpname[MAX_LEN_FILENAME];

int main(int argc, char **argv){  
  int write_phew_particles(char *snap);
  int write_sph_particles(char *snap);  
  int parse_input(int argc, char **argv, char *base_name);

  /* Load HDF5 snapshot */
  /* Inherited from Romeel's specexsnap code */
  /* Should return 1. a gadget header h; 2. Particle data P */
  printf("We have %d parameters here, master.\n", argc);
  strcpy(base_name, argv[1]);
  parse_input(argc, argv, base_name);
  
  load_hdf5(base_name);

  cosmounits();

  /* Obsolete. */
  /* We used to over-plot SPH particles that used to be PhEW. */
  /* BUT, now PhEW almost never becomes SPH particle any more. */
  if(flag_write_sph == 1)
    write_sph_particles(base_name);

  if(flag_write_phew == 1)
    write_phew_particles(base_name);

  /* if(flag_write_halo == 1) */
  /*   write_halo_particles(base_name, haloid);     */

  return 1;
}

int write_sph_particles(char *snap){
  int i;
  char filename[200];
  FILE *outputfile;
  double MeanWeight, LogRho, LogTemp;
  int SfFlag;
  int icount = 0;
  int nskip = h.npart[0] / NSPH;
#ifdef NOSPHSKIP
  nskip = 1;
#endif    
  
  sprintf(filename, "%s.sph", snap);
  outputfile = fopen(filename, "w");
  fprintf(outputfile, "#idx logrho logT Mass MassZ f_wind SfFlag\n");
  for(i=0; i < h.npart[0]; i++) {
    if(i % nskip == 0){ // Selected for output
      /* Density consistent with the rhot.c */
      icount ++;
      LogRho = P[i].Rho * UNIT_M / (pow(UNIT_L, 3) * unit_Density);
      LogRho = log10(LogRho * gheader.HubbleParam * gheader.HubbleParam / OMEGABARYON); 
      MeanWeight = (1 + 4 * XHE) / (1 + P[i].Ne + XHE);
      LogTemp = P[i].Temp * unit_Temp;
      LogTemp = log10(LogTemp * GAMMA_MINUS1 * PROTONMASS / BOLTZMANN * MeanWeight);
      /* Now we have temperature for the particle. */
      SfFlag = (P[i].Sfr > 0) ? 1 : 0;
    
      fprintf(outputfile,"%d %7.5e %6.4f %7.5e %7.5e %5.3e %1d\n",
	      i, LogRho, LogTemp, P[i].Mass, P[i].metal[0], P[i].WindMass, SfFlag);
    } // Mcloud != 0
  }
  fclose(outputfile);
  fprintf(stdout, "Total number of SPH particles: %d\n", icount);
  return 1;
}

/* int is_particle_in_halo(int i){ */
/*   int k; */
/*   for(k = 0; k < 3; k++) */
/*     if(P[i].Pos[k] < HPos) */
/* } */

int write_phew_particles(char *snap){
  int i;
  char filename[200];
  FILE *outputfile;
  double MeanWeight, LogRho, LogTemp, LogRhoa, LogTempa;
  int icount = 0;
  
  sprintf(filename, "%s.phews", snap);
  outputfile = fopen(filename, "w");
  /* 20201104: Add rhoa and Ta */
  fprintf(outputfile, "#idx rhoc Tc rhoa Ta f_cloud f_wind LastSFTime\n");
  for(i=0; i < h.npart[0]; i++) {
    if(P[i].Mcloud > 0){ // Is or Has been a PhEW
      /* Density consistent with the rhot.c */
      icount ++;
      LogRho = P[i].Rho * UNIT_M / (pow(UNIT_L, 3) * unit_Density);
      LogRho = log10(LogRho * gheader.HubbleParam * gheader.HubbleParam / OMEGABARYON); 
      /* Convert internal energy to temperature. From hdf52tipsy. */
      /* From GIZMO:cooling/cooling.c */
      /* static double mhboltz = PROTONMASS / BOLTZMANN; */      
      /* static double yhelium = (1 - HYDROGEN_MASSFRAC) / (4 * HYDROGEN_MASSFRAC); */
      /* T = log10(GAMMA_MINUS1 * u * mhboltz * (1 + 4 * yhelium) / (1 + ne + yhelium)); */
      MeanWeight = (1 + 4 * XHE) / (1 + P[i].Ne + XHE);
      LogTemp = P[i].Temp * unit_Temp;
      LogTemp = log10(LogTemp * GAMMA_MINUS1 * PROTONMASS / BOLTZMANN * MeanWeight);
      /* Now we have temperature for the particle. */
      /* For wind particles, use SFR and Tmax for ambient_density and u_ambient */
      LogRhoa = P[i].Sfr * UNIT_M / (pow(UNIT_L, 3) * unit_Density);
      LogRhoa = log10(LogRhoa * gheader.HubbleParam * gheader.HubbleParam / OMEGABARYON);
      MeanWeight = (1 + 4 * XHE) / (1 + 1.15 + XHE); /* Assuming fully ionized */
      LogTempa = P[i].Tmax;
      LogTempa = log10(LogTempa * GAMMA_MINUS1 * PROTONMASS / BOLTZMANN * MeanWeight);
    
      fprintf(outputfile,"%d %7.5e %6.4f %7.5e %6.4f %5.3f %5.3e %7.5f\n",
	      i, LogRho, LogTemp, LogRhoa, LogTempa, P[i].Mcloud, P[i].WindMass, P[i].LastSFTime);
    } // Mcloud != 0
  }
  fclose(outputfile);
  fprintf(stdout, "Total number of PhEW particles: %d\n", icount);
  return 1;
}

int parse_input(int argc, char **argv, char *base_name)
{
  int i=1;
  while(i < argc) {
    if (!strcmp(argv[i], "-sph"))
      {flag_write_sph=1; ++i;}
    else if (!strcmp(argv[i], "-phew"))
      {flag_write_phew=1; ++i;}
    else if (!strcmp(argv[i], "-halo")){
      flag_write_halo=1;
      ++i;
      haloid = atoi(argv[i]);
      ++i;
    }
    else{
      printf("%s\n", argv[i]);
      strcpy(base_name, argv[i]); ++i;
    }
  }
  return 0;
}
