/* SH190321: Find/Write Particles of a single SKID galaxy */
/* - #idx PID Mass Tmax Age(-1 for gas) */

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
/* #define UNIT_M 3469581.88 */
char grpname[MAX_LEN_FILENAME];
char outname[MAX_LEN_FILENAME];

void main(int argc, char **argv)
{
  FILE *fgalslst;
  int galparts(int gid);
  int i, gid;
  int snapnum;
  char snapstr[4];

  strcpy(basename, "/scratch/shuiyao/data/p50n288sw");
  strcpy(skidbasename, "/scratch/shuiyao/data/p50n288sw");
  snapnum = 108;
  get_snap_string(snapnum, snapstr);
  sprintf(binname, "%s/snap_p50n288gw_%s.bin", basename, snapstr);
  sprintf(auxname, "%s/snap_p50n288gw_%s.aux", basename, snapstr);      
  sprintf(idnumname, "%s/snap_p50n288gw_%s.idnum", basename, snapstr);
  read_tipsy_header(binname);
  read_tipsy_idnum(idnumname);
  read_tipsy_binary(binname, 4);    
  read_tipsy_aux(auxname, 4);  

  if(fgalslst = fopen(argv[1], "r")){ /* First parameter is a file. */
    fprintf(stdout, "Reading Galaxy List File: %s\n", argv[1]);
    i = 0;
    while(1){
      fscanf(fgalslst, "%d\n", &gid);
      fprintf(stdout, "GID = %d\n", gid);
      galparts(gid);
      if(feof(fgalslst)) break;
    }
    fclose(fgalslst);
  }
  else{ /* First parameter is (likely) a number (GID). */
    gid = atoi(argv[1]);
    galparts(gid);
    fprintf(stdout, "Single Galaxy GID = %d\n", gid);
  }
  freeall();
}

int galparts(int gid)
{
  int i, pidx, snapnum;
  FILE *fout;
  int ioffset_star;
  char snapstr[4];  
  /* skidbasename and snapstr are global variable and should have been defined earlier. */
  snapnum = 108;
  get_snap_string(snapnum, snapstr);  
  sprintf(grpname, "%s/gal_z%s.grp", skidbasename, snapstr);
#ifdef GALPARTS_INCLUDE_GAS
  sprintf(outname, "./gal_%d.parts", gid);
#else
  sprintf(outname, "./gal_%d.stars", gid);
#endif  
  read_grp(grpname, 4);
  
  ioffset_star = theader.nsph + theader.ndark;
  fout = fopen(outname, "w");
  fprintf(fout, "#Idx ID Mass Tmax Age\n");
#ifdef GALPARTS_INCLUDE_GAS
  for(i=0;i<theader.nsph;i++){
    if(galid[i] == gid){ /* THAT Gal */
      fprintf(fout, "%d %d %7.5e %5.3f %7.5e\n",
	      i, pids[i],
	      gp[i].mass, aux_gp[i].tmax, -1.0);
    }
  }
#endif
  for(i=0;i<theader.nstar;i++){
    pidx = i + ioffset_star;
    if(galid[pidx] == gid){ /* THAT Gal */
      fprintf(fout, "%8d %8d %7.5e %5.3f %7.5e\n",
	      pidx, pids[pidx],
	      sp[i].mass, aux_sp[i].tmax, aux_sp[i].age);
    }
  }
  
  /* fprintf(stdout, "%s\n", idnumname); */
  /* fprintf(stdout, "%d %d %d\n", pids[0], pids[1], pids[2]); */
  fclose(fout);
}
