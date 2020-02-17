/* SH190321: INCOMPLETE */
/* The original idea is to map the idx to each ID. But the hard thing is that some particles have multiple IDs so it's not straightforward to find the data structure to store the indices. */
  /* - Unlike idx -> ID, where each idx is unique so that an array of IDs is sufficient. With multiple IDs, a single array is not enough. */

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
struct tipsy_header soheader, galheader;
struct star_particle *gals, *sohalos;
char gtpname[MAX_LEN_FILENAME];
char sogtpname[MAX_LEN_FILENAME];
char fofname[MAX_LEN_FILENAME];
char outname[MAX_LEN_FILENAME];

void main()
{
  fprintf(stdout, "UNDER CONSTRUCTION.\n");
}

int skidsopid()
{
  int i, j, k, pass;
  double dx, dr;
  double mstar, rmin, hsml;
  int snapnum, numngb, nn;
  char snapstr[4];
  FILE *fout;
  /* 1. r within a few r_srh */
  /* 2. Mstar within 1.0 dex */

  strcpy(basename, "/scratch/shuiyao/data/p50n288zw");
  strcpy(skidbasename, "/scratch/shuiyao/data/p50n288zw");
  snapnum = 108;
  get_snap_string(snapnum, snapstr);
  sprintf(gtpname, "%s/gal_z%s.gtp", skidbasename, snapstr);
  sprintf(sogtpname, "%s/so_z%s.sogtp", skidbasename, snapstr);
  sprintf(binname, "%s/snap_p50n288gw_%s.bin", basename, snapstr);    
  sprintf(idnumname, "%s/snap_p50n288gw_%s.idnum", basename, snapstr);
  sprintf(outname, "%s/pids_z%s.grp", basename, snapstr);    
  /* read_gtp(gtpname); */
  /* read_sogtp(sogtpname); */
  read_tipsy_header(binname);
  read_tipsy_idnum(idnumname, 0);

  fout = fopen(outname, "w");
  fprintf(fout, "#idx nn mfof mvir mstar rmin numngb\n");
  
  /* fprintf(stdout, "%s\n", idnumname);   */
  /* fprintf(stdout, "%d %d %d\n", pids[0], pids[1], pids[2]); */
  freeall();
  /* fclose(fout); */
}
