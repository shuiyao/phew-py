/* SH191129: Select a certain halo from a snapshot */
/* Usage: get_halo_particles snapname sogrpname haloid */

/* - #idx PID Mass Pos[0] Pos[1] Pos[2] Vel[0] Vel[1] Vel[2] Density Temperature Z */

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include "tipsydefs.h"
#include "gadget3io.h"

#define BOXSIZE 1.0
#define HALFBOX 0.5
/* #define UNIT_M 433697.735404 */
#define UNIT_M 3469581.88

/* char grpname[MAX_LEN_FILENAME]; */
char sogrpname[MAX_LEN_FILENAME];
char snapname[MAX_LEN_FILENAME];
char binoutname[MAX_LEN_FILENAME];
char auxoutname[MAX_LEN_FILENAME];
char idnumoutname[MAX_LEN_FILENAME];
int flag_allstars = 0; /* if 1, dump all stars. Otherwise, only dump stars from a certain Mvir range. */
struct tipsy_header theaderout;

void main(int argc, char **argv)
{
  int i, gid, hid;
  int snapnum;

  hid = atoi(argv[3]);
  strcpy(sogrpname, argv[2]);
  strcpy(snapname, argv[1]);
  
  sprintf(binname, "%s.bin", snapname);
  sprintf(auxname, "%s.aux", snapname);
  sprintf(idnumname, "%s.idnum", snapname);
  sprintf(binoutname, "%s_%s.bin", snapname, argv[3]);
  sprintf(auxoutname, "%s_%s.aux", snapname, argv[3]);
  sprintf(idnumoutname, "%s_%s.idnum", snapname, argv[3]);

  fprintf(stdout, "Binnary File: %s\n", binname); /* File check. */
  
  read_tipsy_header(binname);
  read_tipsy_idnum(idnumname, 4);
  read_tipsy_binary(binname, 4);
  read_tipsy_aux(auxname, 4);
  /* read_grp(grpname, 4); */
  read_sogrp(sogrpname, 4);

  int pidx;
  int Nsph, Ndark, Nstar, Noffset;
  FILE *fbinout, *fauxout, *fidnumout;
  fbinout = fopen(binoutname, "w");
  fauxout = fopen(auxoutname, "w");
  fidnumout = fopen(idnumoutname, "w");
  
  /* How many gas particles are in the hid halo? */
  Nsph = 0;
  Ndark = 0;
  Nstar = 0;
  Noffset = 0;
  for(i=0;i<theader.nsph;i++)
    if(sohid[Noffset+i] == hid){
      Nsph ++; 
    }
  Noffset = theader.nsph;
  for(i=0;i<theader.ndark;i++)
    if(sohid[Noffset+i] == hid){
      Ndark ++; 
    }
  Noffset = theader.nsph + theader.ndark;
  for(i=0;i<theader.nstar;i++)
    if(sohid[Noffset+i] == hid){
      Nstar ++; 
    }
  theaderout.time = theader.time;
  theaderout.nbodies = Nsph + Ndark + Nstar;
  theaderout.ndim = 3;
  theaderout.nsph = Nsph;
  theaderout.ndark = Ndark;
  theaderout.nstar = Nstar;

  fprintf(stdout, "Halo ID: %d\n", hid);
  fprintf(stdout, "--------------------------------\n");
  fprintf(stdout, "Nsph = %d\n", Nsph);
  fprintf(stdout, "Ndark = %d\n", Ndark);
  fprintf(stdout, "Nstar = %d\n", Nstar);    

  fwrite(&theaderout, sizeof(theaderout), 1, fbinout);
  /* Read and Write */
  Noffset = 0;
  for(i=0;i<theader.nsph;i++)
    if(sohid[Noffset+i] == hid){
      fwrite(&gp[i], sizeof(struct gas_particle), 1, fbinout);
      fwrite(&aux_gp[i], sizeof(struct aux_gas_data), 1, fauxout);
      fwrite(&pids[Noffset+i], sizeof(int), 1, fidnumout);      
    }
  Noffset = theader.nsph;
  for(i=0;i<theader.ndark;i++)
    if(sohid[Noffset+i] == hid){
      fwrite(&dp[i], sizeof(struct dark_particle), 1, fbinout);
      fwrite(&pids[Noffset+i], sizeof(int), 1, fidnumout);      
    }
  Noffset = theader.nsph + theader.ndark;
  for(i=0;i<theader.nstar;i++)
    if(sohid[Noffset+i] == hid){
      fwrite(&sp[i], sizeof(struct star_particle), 1, fbinout);
      fwrite(&aux_sp[i], sizeof(struct aux_star_data), 1, fauxout);
      fwrite(&pids[Noffset+i], sizeof(int), 1, fidnumout);      
    }
  
  fclose(fbinout);
  fclose(fauxout);
  fclose(fidnumout);
  freeall();
}

