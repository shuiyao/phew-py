#ifndef _GADGET3IO_H
#define _GADGET3IO_H

#include <stdio.h>

#define MAX_LEN_FILENAME 500
#define MAX_NSKID_GRPS 50000

/* Variables */
extern char basename[MAX_LEN_FILENAME];
extern char skidbasename[MAX_LEN_FILENAME];
extern char binname[MAX_LEN_FILENAME], auxname[MAX_LEN_FILENAME], idnumname[MAX_LEN_FILENAME];
extern char grpname[MAX_LEN_FILENAME], gtpname[MAX_LEN_FILENAME];
extern char sogrpname[MAX_LEN_FILENAME], sogtpname[MAX_LEN_FILENAME];
extern char fofname[MAX_LEN_FILENAME];

extern int *pids;
extern int *pidxs;
extern int *sohid;
extern int *galid;

extern struct gas_particle *gp;
extern struct dark_particle *dp;
extern struct star_particle *sp;
extern struct aux_gas_data *aux_gp;
extern struct aux_star_data *aux_sp;
extern struct star_particle *gals, *sohalos;
extern struct tipsy_header theader, soheader, galheader;
extern struct fof_group{
  int gid;
  float mstar;
  float mgal;
  float mfof;
  double pos[3];
}
  *fofgrps;
extern struct skid_group {
  int gid;
  float mgas;
  float mstar;
  double pos[3];
  double metals[4];
  double sfr;
}
  *skid_grps;

extern int Nfofgrps;

/* Functions */
void get_snap_string(int snapnum, char *snapstr);
void read_header(struct tipsy_header *header, FILE *ftipsy);
int headerinfo(struct tipsy_header *header);
int read_tipsy_header(char *filename);
int read_tipsy_binary(char *filename, int typecode);
int read_tipsy_aux(char *filename, int typecode);
int read_tipsy_idnum(char *filename, int reverse_map);
int read_gtp(char *filename);
int read_sogtp(char *filename);
int read_sogrp(char *filename, int typecode);
int read_grp(char *filename, int typecode);
int read_fofgroup(char *filename);
void freeall();

#endif
