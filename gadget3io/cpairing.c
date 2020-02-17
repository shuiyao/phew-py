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
#define RSEARCH 0.001
#define MS_TOLERANCE 1.0
struct tipsy_header soheader, galheader;
struct star_particle *gals, *sohalos;
char gtpname[MAX_LEN_FILENAME];
char sogtpname[MAX_LEN_FILENAME];
char fofname[MAX_LEN_FILENAME];

void main()
{
  int find_fof_so_pairs();
  find_fof_so_pairs();
}

int find_fof_so_pairs()
{
  int i, j, k, pass;
  double dx, dr;
  double mstar, rmin, hsml;
  int snapnum, numngb, nn;
  char snapstr[4];
  FILE *fout;
  /* 1. r within a few r_srh */
  /* 2. Mstar within 1.0 dex */
  fout = fopen("fofpairs.dat", "w");
  fprintf(fout, "#idx nn mfof mvir mstar rmin numngb\n");

  strcpy(basename, "/scratch/shuiyao/data/p50n288zw");
  strcpy(skidbasename, "/scratch/shuiyao/data/p50n288zw");
  snapnum = 108;
  get_snap_string(snapnum, snapstr);
  sprintf(gtpname, "%s/gal_z%s.gtp", skidbasename, snapstr);
  sprintf(sogtpname, "%s/so_z%s.sogtp", skidbasename, snapstr);
  sprintf(fofname, "%s/FoFGroups_%s/fofgroups.dat", skidbasename, snapstr);
  read_gtp(gtpname);
  read_sogtp(sogtpname);
  Nfofgrps = read_fofgroup(fofname);
  
  for(i=0;i<Nfofgrps;i++){
    numngb = 0; hsml = 0; nn = -1; rmin = BOXSIZE;
    if((float)(i) / 1000.0 == (i / 1000))
      fprintf(stdout, "(%5d/%5d) FoF groups done\n", i, Nfofgrps);
    while(numngb < 1 && hsml <= 10.0*RSEARCH){
      hsml += RSEARCH;
      for(j=0;j<soheader.nstar;j++){
	dr = 0; pass = 1;
	for(k=0;k<3;k++){
	  dx = fofgrps[i].pos[k] - sohalos[j].pos[k];
	  if(dx > HALFBOX) dx = BOXSIZE - dx;
	  if(dx > hsml){pass = 0; break;}
	  dr += dx * dx;
	}
	if(pass == 0) continue;
	dr = sqrt(dr);
	if(dr > hsml) continue;
	mstar = log10(gals[j].mass * UNIT_M * 1.e10 / 0.7);
	if(fabs(fofgrps[i].mstar - mstar) > MS_TOLERANCE) continue;
	// Successful Match.
	numngb += 1;
	if(rmin > dr) {rmin = dr; nn = j;}
      } // Sohalos
    } // while
    if(numngb > 0)
      fprintf(fout, "%5d %5d %6.3f %6.3f %6.3f %7.5f %2d\n",
	      i, nn, fofgrps[i].mfof,
	      log10(sohalos[nn].mass * UNIT_M * 1.e10 / 0.7),
	      log10(gals[nn].mass * UNIT_M * 1.e10 / 0.7),
	      rmin, numngb);
    else
      fprintf(fout, "%5d %5d %5.3f %5.3f %5.3f %7.5f %2d\n",
	      i, nn, fofgrps[i].mfof,
	      0.0,
	      0.0,
	      rmin, numngb);
  } // Nfofgrps
  freeall();
  fclose(fout);
}
