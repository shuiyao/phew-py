/* 180630 [CRITICAL] Bug Fix. Fix the mass and other properties for stars. */

/* Need Mass Msub, Rsub, WindMass, dR, Mcloud */

#include <stdio.h>
#include "gadgetdefs.h"
#include "tipsydefs.h"
#include "loadhdf5.h"
#include <hdf5.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define NSPH 50000
#define BOXSIZE 25000
#define UNIT_MASS 433697.735404
#define OMEGABARYON 0.045
#define HUBBLEPARAM 0.7

char modelname[MAX_LEN_FILENAME];
char base_name[MAX_LEN_FILENAME], skid_base_name[MAX_LEN_FILENAME];
char binname[MAX_LEN_FILENAME], hdf5name[MAX_LEN_FILENAME];
char grpname[MAX_LEN_FILENAME], gtpname[MAX_LEN_FILENAME];
char sogrpname[MAX_LEN_FILENAME], sogtpname[MAX_LEN_FILENAME];
char soparname[MAX_LEN_FILENAME];

int snapnum;
char snapstr[4];

int main(int argc, char **argv){  
  /* Load HDF5 snapshot */
  /* Inherited from Romeel's specexsnap code */
  /* Should return 1. a gadget header h; 2. Particle data P */
  void get_halo_particles();
  
  printf("We have %d parameters here, master.\n", argc);
  strcpy(modelname, argv[1]);
  snapnum = atoi(argv[2]);

  strcat(strcpy(base_name, "/proj/shuiyao/"), modelname);
  strcat(strcpy(skid_base_name, "/proj/shuiyao/"), modelname);
  get_snap_string(snapnum, snapstr);  
  
  strcpy(skid_base_name, base_name);
  sprintf(hdf5name, "%s/snapshot_%s", base_name, snapstr);
  sprintf(binname, "%s/snapshot_%s.bin", base_name, snapstr);
  fprintf(stdout, "HDF5 File: %s\n", hdf5name); /* File check. */
  sprintf(grpname, "%s/gal_z%s.grp", skid_base_name, snapstr);
  sprintf(sogrpname, "%s/so_z%s.sogrp", skid_base_name, snapstr);
  sprintf(soparname, "%s/so_z%s.par", skid_base_name, snapstr);  

  load_hdf5(hdf5name);

  cosmounits();

  read_tipsy_header(binname); // Get theader;
  read_grp(grpname, 0); // Get galid
  read_sogrp(sogrpname, 0);  // Get sohid
  read_sopar(soparname);   
  get_halo_particles();

  free(hosthid);  
  free(sohid);
  free(galid);
  free(P);

  return 1;
}

void get_halo_particles(){
  int i, j, hid;
  FILE *foutm11, *foutm12, *foutm13;
  char outm11[MAX_LEN_FILENAME];
  char outm12[MAX_LEN_FILENAME];
  char outm13[MAX_LEN_FILENAME];
  char printline[200];
  double mvir, rvir;
  double dx, dr;

  sprintf(sogtpname, "%s/so_z%s.sogtp", skid_base_name, snapstr);
  read_sogtp(sogtpname);
  for(i=0;i<soheader.nstar;i++){ /* Convert Msub to log(msolar) */
    if(sohalos[i].mass > 0) sohalos[i].mass = log10(sohalos[i].mass * UNIT_MASS / HUBBLEPARAM) + 10.;
  }
  sprintf(outm11, "/scratch/shuiyao/scidata/gadget3io/%s/%s_%s.gas.mh11", modelname, modelname, snapstr);
  sprintf(outm12, "/scratch/shuiyao/scidata/gadget3io/%s/%s_%s.gas.mh12", modelname, modelname, snapstr);
  sprintf(outm13, "/scratch/shuiyao/scidata/gadget3io/%s/%s_%s.gas.mh13", modelname, modelname, snapstr);    
  foutm11 = fopen(outm11, "w");
  foutm12 = fopen(outm12, "w");
  foutm13 = fopen(outm13, "w");
  fprintf(foutm11, "#Idx GID HID Mass WMass Mc dr Rvir\n");
  fprintf(foutm12, "#Idx GID HID Mass WMass Mc dr Rvir\n");
  fprintf(foutm13, "#Idx GID HID Mass WMass Mc dr Rvir\n");
  for(i=0;i<theader.nsph;i++){
    hid = sohid[i];
    if(hid <= 0) continue;
    if(hosthid[hid] != hid) continue; // Not a central halo
    if(hid > 0){
      mvir = sohalos[hid-1].mass;
      rvir = sohalos[hid-1].eps / HUBBLEPARAM * BOXSIZE * theader.time;
      dr = 0.0;
      for(j=0;j<3;j++){
	dx = (sohalos[hid-1].pos[j] + 0.5) * BOXSIZE - P[i].Pos[j];
	dx = (dx > BOXSIZE / 2) ? BOXSIZE - dx : dx;
	dr += dx * dx;
      }
      dr = sqrt(dr) * theader.time / HUBBLEPARAM;
#ifdef PHEW_EXTRA_OUTPUT      
      sprintf(printline, "%8d %5d %5d %7.5e %7.5e %5.3f %7.2f %7.3f\n",
	      i, galid[i], sohid[i],
	      P[i].Mass, P[i].WindMass,
	      P[i].Mcloud, dr, rvir
	      );
#else
      sprintf(printline, "%8d %5d %5d %7.5e %7.5e %5.3f %7.2f %7.3f\n",
	      i, galid[i], sohid[i],
	      P[i].Mass, 0.0,
	      P[i].DelayTime, dr, rvir
	      );
#endif      
      /* sp[i].mass, fabs(aux_sp[i].tmax), aux_sp[i].age); */
      if((11.0 < mvir) && (mvir < 11.5)) fprintf(foutm11, printline);
      else if((11.85 < mvir) && (mvir < 12.15)) fprintf(foutm12, printline);
      else if((12.85 < mvir) && (mvir < 13.15)) fprintf(foutm13, printline);	
    }
  }
  fclose(foutm11);
  fclose(foutm12);
  fclose(foutm13);
}

