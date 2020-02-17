/* SH190423: Compile information for all stars from a snapshot */
/* Based on the galparts.c */
/* Usage: allstars modelname snapnum all */
/* - #idx PID GID HID Mass Tmax Age */

/* **************** WARNING on the Tmax **************** */
/* In versions after p50n288dsw, Tmax is negative if the last time a particle is SF is when it is launched in a wind. But to be consistent, I ignore this information now and enforce all Tmax to be positive by having fabs(auxsp.tmax). */

/* I changed it back to include NEGATIVE tmax for p50n288fiducial */

/* **************** ******************* **************** */

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include "tipsydefs.h"
#include "/scratch/shuiyao/sci/gadget3io/gadget3io.h"

#define BOXSIZE 1.0
#define HALFBOX 0.5
/* #define UNIT_M 433697.735404 */
#define UNIT_M 3469581.88

/* char grpname[MAX_LEN_FILENAME]; */
/* char sogrpname[MAX_LEN_FILENAME]; */
int flag_allstars = 0; /* if 1, dump all stars. Otherwise, only dump stars from a certain Mvir range. */

void main(int argc, char **argv)
{
  FILE *fgalslst;
  int i, gid, hid;
  int snapnum;
  char snapstr[4];
  char modelname[20];

  snapnum = atoi(argv[2]);
  /* modelname is argv[1] */
  strcpy(modelname, argv[1]);
  if((argc == 4) && (!strcmp(argv[3], "all"))){
    fprintf(stdout, "Doing every single star.\n");
    flag_allstars = 1;
  }
  
  strcat(strcpy(basename, "/scratch/shuiyao/data/"), modelname);
  strcat(strcpy(skidbasename, "/scratch/shuiyao/data/"), modelname);
  get_snap_string(snapnum, snapstr);
  sprintf(binname, "%s/snap_p50n288gw_%s.bin", basename, snapstr);
  sprintf(auxname, "%s/snap_p50n288gw_%s.aux", basename, snapstr);
  sprintf(idnumname, "%s/snap_p50n288gw_%s.idnum", basename, snapstr);
  fprintf(stdout, "Binnary File: %s\n", binname); /* File check. */
  sprintf(grpname, "%s/gal_z%s.grp", skidbasename, snapstr);
  sprintf(sogrpname, "%s/so_z%s.sogrp", skidbasename, snapstr);
  
  read_tipsy_header(binname);
  read_tipsy_idnum(idnumname, 0);
  read_tipsy_binary(binname, 4);
  read_tipsy_aux(auxname, 4);
  read_grp(grpname, 4);
  read_sogrp(sogrpname, 4);

  int ioffset_star;
  FILE *fout;
  int pidx;
  char outname[MAX_LEN_FILENAME];
  FILE *foutm11, *foutm12, *foutm13;
  char outm11[MAX_LEN_FILENAME];
  char outm12[MAX_LEN_FILENAME];
  char outm13[MAX_LEN_FILENAME];
  char printline[200];
  double mvir;

  if(flag_allstars){
    ioffset_star = theader.nsph + theader.ndark;
    sprintf(outname, "/scratch/shuiyao/scidata/gadget3io/%s/%s_%s.stars", modelname, modelname, snapstr);
    fprintf(stdout, "Writing file: %s\n", outname);
    fout = fopen(outname, "w");
    fprintf(fout, "#Idx ID GID HID Mass Tmax Age\n");
    for(i=0;i<theader.nstar;i++){
      pidx = i + ioffset_star;
      fprintf(fout, "%8d %8d %5d %5d %7.5e %5.3f %7.5e\n",
	      pidx, pids[pidx], galid[pidx], sohid[pidx],
	      sp[i].mass, aux_sp[i].tmax, aux_sp[i].age);	      
	      /* sp[i].mass, fabs(aux_sp[i].tmax), aux_sp[i].age); */
    }
    fclose(fout);
  }
  else{ /* flag_allstars = False */
    ioffset_star = theader.nsph + theader.ndark;
    /* fprintf(stdout, "%s\n", sogtpname);   */
    sprintf(sogtpname, "%s/so_z%s.sogtp", skidbasename, snapstr);
    read_sogtp(sogtpname);
    for(i=0;i<soheader.nstar;i++){ /* Convert Msub to log(msolar) */
      if(sohalos[i].mass > 0) sohalos[i].mass = log10(sohalos[i].mass * UNIT_M / 0.7) + 10.;
    }
    sprintf(outm11, "/scratch/shuiyao/scidata/gadget3io/%s/%s_%s.stars.mh11", modelname, modelname, snapstr);
    sprintf(outm12, "/scratch/shuiyao/scidata/gadget3io/%s/%s_%s.stars.mh12", modelname, modelname, snapstr);
    sprintf(outm13, "/scratch/shuiyao/scidata/gadget3io/%s/%s_%s.stars.mh13", modelname, modelname, snapstr);    
    foutm11 = fopen(outm11, "w");
    foutm12 = fopen(outm12, "w");
    foutm13 = fopen(outm13, "w");
    fprintf(foutm11, "#Idx ID GID HID Mass Tmax Age\n");
    fprintf(foutm12, "#Idx ID GID HID Mass Tmax Age\n");
    fprintf(foutm13, "#Idx ID GID HID Mass Tmax Age\n");
    for(i=0;i<theader.nstar;i++){
      pidx = i + ioffset_star;
      hid = sohid[pidx];
      if(hid > 0){
	mvir = sohalos[hid-1].mass;
	sprintf(printline, "%8d %8d %5d %5d %7.5e %5.3f %7.5e\n",
		pidx, pids[pidx], galid[pidx], sohid[pidx], 
		sp[i].mass, fabs(aux_sp[i].tmax), aux_sp[i].age);
	if((11.0 < mvir) && (mvir < 11.5)) fprintf(foutm11, printline);
	else if((11.85 < mvir) && (mvir < 12.15)) fprintf(foutm12, printline);
	else if((12.85 < mvir) && (mvir < 13.15)) fprintf(foutm13, printline);	
      }
    }
    fclose(foutm11);
    fclose(foutm12);
    fclose(foutm13);    
  }
  freeall();
}

