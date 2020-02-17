#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include "tipsydefs.h"
#include "gadget3io.h"

#ifdef MARK64
void read_header(struct tipsy_header *, FILE *);
#endif

char basename[MAX_LEN_FILENAME], skidbasename[MAX_LEN_FILENAME];
char binname[MAX_LEN_FILENAME], auxname[MAX_LEN_FILENAME], idnumname[MAX_LEN_FILENAME];
char grpname[MAX_LEN_FILENAME], gtpname[MAX_LEN_FILENAME];
char sogrpname[MAX_LEN_FILENAME], sogtpname[MAX_LEN_FILENAME];
char fofname[MAX_LEN_FILENAME];
struct gas_particle *gp;
struct dark_particle *dp;
struct star_particle *sp;
struct aux_gas_data *aux_gp;
struct aux_star_data *aux_sp;
struct tipsy_header theader, soheader, galheader;
struct star_particle *gals, *sohalos;
struct skid_group *skid_grps;
struct fof_group *fofgrps;
int *pids; /* Map from pidx -> PID */
int *pidxs; /* Map from PID -> pidx */
int *sohid; /* Map from pidx -> HID */
int *galid; /* Map from pidx -> GID */

int Nfofgrps;

int test(int argc, char **argv)
{
  int snapnum;
  int i, j, ngals;
  char snapstr[4];
  
  strcpy(basename, "./data");
  strcpy(skidbasename, "./data");
  snapnum = 36;

  get_snap_string(snapnum, snapstr);
  sprintf(binname, "%s/snap_p25n144phew_%s.bin", basename, snapstr);
  sprintf(auxname, "%s/snap_p25n144phew_%s.aux", basename, snapstr);
  sprintf(grpname, "%s/gal_z%s.grp", skidbasename, snapstr);
  sprintf(gtpname, "%s/gal_z%s.gtp", skidbasename, snapstr);
  sprintf(sogrpname, "%s/so_z%s.sogrp", skidbasename, snapstr);
  sprintf(sogtpname, "%s/so_z%s.sogtp", skidbasename, snapstr);
  /* sprintf(fofname, "%s/FoFGroups_%s/fofgroups.dat", skidbasename, snapstr); */
  sprintf(fofname, "%s/FoFGroups_058/fofgroups.dat", skidbasename, snapstr);

  read_tipsy_binary(binname, 4);
  read_tipsy_aux(auxname, 4);
  read_gtp(gtpname);
  read_sogtp(sogtpname);
  /* ngals = read_skid(snapnum); */
  Nfofgrps = read_fofgroup(fofname);
  freeall();
  return 1;
}

void freeall()
{
  free(gp);
  free(dp);
  free(sp);
  free(aux_gp);
  free(aux_sp);
  free(gals);
  free(sohalos);
  free(skid_grps);
  free(pids);  
  /* free(pidxs); */
}

int read_tipsy_header(char *filename){
  FILE *fp;
  if(!(fp = fopen(filename, "r")))
    {printf("Can't open file %s! Force exit.\n", filename); exit(-1);}    
#ifdef MARK64
  printf("I am reading for 64 machine!\n");
  read_header(&theader, fp);
#else
  printf("I am reading for 32 machine!\n");
  if(fread((char *)&theader,sizeof(theader),1,fp) == 0)
    {printf("BREAK!\n"); exit(-1);}
#endif
  //  headerinfo(&header);
  printf("*** Header Information ***\n");
  printf("header.time: %f\n", theader.time);
  printf("header.nbodies: %d\n", theader.nbodies);
  printf("header.ndim: %d\n", theader.ndim);
  printf("header.nsph: %d\n", theader.nsph);
  printf("header.ndark: %d\n", theader.ndark);
  printf("header.nstar: %d\n", theader.nstar);
  printf("**************************\n");
  fclose(fp);
}

int read_tipsy_binary(char *filename, int typecode)
{
  FILE *fp;
  if(!(fp = fopen(filename, "r")))
    {printf("Can't open file %s! Force exit.\n", filename); exit(-1);}    
#ifdef MARK64
  printf("I am reading for 64 machine!\n");
  read_header(&theader, fp);
#else
  printf("I am reading for 32 machine!\n");
  if(fread((char *)&theader,sizeof(theader),1,fp) == 0)
    {printf("BREAK!\n"); exit(-1);}
#endif
  //  headerinfo(&header);
  printf("*** Header Information ***\n");
  printf("header.time: %f\n", theader.time);
  printf("header.nbodies: %d\n", theader.nbodies);
  printf("header.ndim: %d\n", theader.ndim);
  printf("header.nsph: %d\n", theader.nsph);
  printf("header.ndark: %d\n", theader.ndark);
  printf("header.nstar: %d\n", theader.nstar);
  printf("**************************\n");
  if(theader.nsph != 0) {
    gp = (struct gas_particle *)
      malloc(theader.nsph*sizeof(*gp));
    if(gp == NULL) {
      printf("<sorry, no memory for GAS particles, master>\n") ;
      return -1;
    }
  }
  if(theader.ndark != 0 && typecode >= 1) {
    dp = (struct dark_particle *)
      malloc(theader.ndark*sizeof(*dp));
    if(dp == NULL) {
      printf("<sorry, no memory for DARK particles, master>\n") ;
      return -1;
    }
  }
  if(theader.nstar != 0 && typecode == 4) {
    sp = (struct star_particle *)
      malloc(theader.nstar*sizeof(*sp));
    if(sp == NULL) {
      printf("<sorry, no memory for STAR particles, master>\n") ;
      return -1;
    }
  }
  fread((char *)gp,sizeof(struct gas_particle),	theader.nsph,fp) ;
  if(typecode >= 1)
    fread((char *)dp,sizeof(struct dark_particle), theader.ndark,fp) ;
  if(typecode == 4)
    fread((char *)sp,sizeof(struct star_particle), theader.nstar,fp) ;
  fclose(fp);
  printf("---> Reading .bin done.\n");
  fprintf(stdout,"read time %f\n",theader.time);
  return 1;
}

int read_tipsy_aux(char *filename, int typecode)
{
  FILE *fp;
  if(!(fp = fopen(filename, "r")))
    {printf("Can't open file %s! Force exit.\n", filename); exit(-1);}
/* Read in auxiliary particle data file */
  aux_gp = (struct aux_gas_data *) malloc(theader.nsph*sizeof(*aux_gp));
  if(aux_gp == NULL) {
    printf("<sorry, no memory for AUX_GP particles, master>\n") ;
    return -1;
  }
  fread(aux_gp,sizeof(struct aux_gas_data),theader.nsph,fp);
  if(typecode == 4){
    aux_sp = (struct aux_star_data *) malloc(theader.nstar*sizeof(*aux_sp));
    if(aux_sp == NULL) {
      printf("<sorry, no memory for AUX_SP particles, master>\n") ;
      return -1;
    }
    fread(aux_sp,sizeof(*aux_sp),theader.nstar,fp);
  }
  fclose(fp);
  fprintf(stdout, "---> Reading .aux done.\n");  
  return 1;
}

/* Create the pidx -> PID map */
/* Turn on the reverse_map will also create the PID -> pidx map */
/*   - Assume that theader.nbodies and int *pids are already established */
/*   - Note: Each pidx should have a unique PID, though some star PID is in the SPH regime */
int read_tipsy_idnum(char *filename, int reverse_map)
{
  int i;
  unsigned char *buf;
  FILE *fp;

  if(!(fp = fopen(filename, "r")))
    {printf("Can't open file %s! Force exit.\n", filename); exit(-1);}    
  
  pids = (int *) malloc(theader.nbodies*sizeof(int));
  if(pids == NULL) {
    printf("<sorry, no memory for PID, master>\n") ;
    return -1;
  }
  fread(pids,sizeof(int),theader.nbodies,fp);    
  // ID is always one line more!!
  fclose(fp);
  fprintf(stdout, "read_tipsy_idunm(): Number of particles: %ld\n", theader.nbodies);
  if(reverse_map == 1){
    pidxs = (int *) malloc((theader.nbodies+1)*sizeof(int));
    pidxs[0] = -1; // dummy, no PID == 0
    for(i=0; i<theader.nbodies; i++)
      pidxs[pids[i]] = i;
    }
  printf("---> Reading .idnum done.\n");
}

int read_gtp(char *filename)
{
  FILE *fp;
  fp = fopen(filename, "r");
#ifdef MARK64
  printf("I am reading for 64 machine!\n");
  read_header(&galheader, fp);
#else
  printf("I am reading for 32 machine!\n");
  if(fread((char *)&galheader,sizeof(struct tipsy_header),1,fp) == 0) // PAD HERE
    {printf("EMPTY SOGTP FILE, BREAK!\n"); exit(-1);}
#endif
/*   headerinfo(&header); */
  printf("*** gtp Information ***\n");
  printf("galheader.nstar: %d\n", galheader.nstar);
  printf("**************************\n");

  if(galheader.nsph != 0) {fprintf(stdout, "BAD: NSPH != 0\n"); exit(-1);}
  if(galheader.ndark != 0) {fprintf(stdout, "BAD: NDARK != 0\n"); exit(-1);}
  if(galheader.nstar != 0) {
    gals = (struct star_particle *)
      malloc(galheader.nstar*sizeof(*gals));
    if(gals == NULL) {
      printf("<sorry, no memory for GALS particles, master>\n") ;
      return -1;
    }
    fread((char *)gals,sizeof(struct star_particle),
	  galheader.nstar,fp) ;
  }
  else{
    fprintf(stdout, "Warning: No galaxy in this snapshot!\n");
  }
  fclose(fp);
  printf("---> Reading .gtp done.\n");
  return 1;
}

int read_sogtp(char *filename)
{
  FILE *fp;
  fp = fopen(filename, "r");
#ifdef MARK64
  printf("I am reading for 64 machine!\n");
  read_header(&soheader, fp);
#else
  printf("I am reading for 32 machine!\n");
  if(fread((char *)&soheader,sizeof(struct tipsy_header)-4,1,fp) == 0) // PAD HERE
    {printf("EMPTY SOGTP FILE, BREAK!\n"); exit(-1);}
#endif
/*   headerinfo(&header); */
  printf("*** Sogtp Information ***\n");
  printf("soheader.nstar: %d\n", soheader.nstar);
  printf("**************************\n");

  if(soheader.nsph != 0) {fprintf(stdout, "BAD: NSPH != 0\n"); exit(-1);}
  if(soheader.ndark != 0) {fprintf(stdout, "BAD: NDARK != 0\n"); exit(-1);}
  if(soheader.nstar != 0) {
    sohalos = (struct star_particle *)
      malloc(soheader.nstar*sizeof(*sohalos));
    if(sohalos == NULL) {
      printf("<sorry, no memory for SOHALOS particles, master>\n") ;
      return -1;
    }
  }
  fread((char *)sohalos,sizeof(struct star_particle),
  	soheader.nstar,fp) ;
  fclose(fp);
  printf("---> Reading .sogtp done.\n");
  return 1;
}

int read_fofgroup(char *filename)
{
  FILE *fp;
  int i, j;
  char dummyline[500];
  int idummy;
  float fdummy;
  float mfof, mstar, mgal, sfr;
  double pos[3];
  int ngrps = 0;
  fofgrps = (struct fof_group *) malloc(MAX_NSKID_GRPS*sizeof(*fofgrps));
  i = 0;
  if(!(fp = fopen(filename, "r"))){
    fprintf(stderr, "Can't open FoF file! Force exit.\n"); exit(-1);}
  fscanf(fp, "%[^\n]", dummyline); // First line, header
  // #CpuID GroupID Mfof Mstar Mgal SFR x y z    
  while(1){
    fscanf(fp, "%d %d %f %f %f %f %lf %lf %lf", 
	   &idummy, &idummy, &mfof, &mstar, &mgal,
	   &sfr, &pos[0], &pos[1], &pos[2]);
    // We ignore dark halos (no stars / no gas)
    if(mgal >= 0){
      fofgrps[ngrps].gid = ngrps;
      fofgrps[ngrps].mgal = mgal;
      fofgrps[ngrps].mstar = mstar;
      fofgrps[ngrps].mfof = mfof;
      for (j=0;j<3;j++) fofgrps[ngrps].pos[j] = pos[j];
      ngrps ++;
      if(ngrps >= MAX_NSKID_GRPS){
	fprintf(stderr, "Number of FoF Groups exceeds limit.\n"); exit(-1);}
    }
    if(feof(fp)) break;
    i++;
  }
  fclose(fp);
  return ngrps;
}

int read_skid(int snapnum)
{
  FILE *inputfile;
  char skidname[200], grpname[200];
  int n, i, j, gid;
  char snapstr[4];
  char dummyline[500];
  double starmass, totmass, gasmass;
  int nskid_grps;

  get_snap_string(snapnum, snapstr);
  strcat(strcpy(skidname, skidbasename), snapstr);
  strcat(skidname, ".stat");
  strcat(strcpy(grpname, skidbasename), snapstr);
  strcat(grpname, ".grp");
  printf("Skidname: %s\n", skidname);
  printf("Grpname: %s\n", grpname);
  // Reading the .stat file
  skid_grps = (struct skid_group *) malloc(MAX_NSKID_GRPS*sizeof(*skid_grps));
  i = 0;
  if(!(inputfile = fopen(skidname, "r")))
    {printf("Can't open .stat file! Force exit.\n"); exit(-1);}
  while(1){
    fscanf(inputfile, "%d %d %f %f %f %[^\n]", 
	   &gid, &n, &totmass, &gasmass, &starmass, dummyline);
    skid_grps[i].gid = gid;
    skid_grps[i].mgas = gasmass;
    skid_grps[i].mstar = starmass;
    skid_grps[i].mgas = 0.;
    skid_grps[i].mstar = 0.;
    skid_grps[i].sfr = 0.;
    for (j=0;j<4;j++) skid_grps[i].metals[j] = 0.0;
    if(feof(inputfile)) break;
    i++;
  }
  nskid_grps = gid; // SHOULD BE THE LARGEST GID
  fclose(inputfile);
  printf("Reading .stat DONE. %d=%d groups read.\n", i, nskid_grps);
  /* // BEGIN READING .GRP FILE */
  /* if(!(inputfile = fopen(grpname, "r"))) */
  /*   {printf("Can't open .grp file! Force exit.\n"); exit(-1);} */
  /* fscanf(inputfile, "%d %[^\n]", &gid, dummyline); // Dummy Line */
  /* printf("First Line in .grp: %d\n", gid); */
  /* for(i=0;i<theader.nsph;i++){ // Add up Gas Metals */
  /*   fscanf(inputfile, "%d\n", &gid); */
  /* } */
  /* for(i=0;i<theader.ndark;i++){ // Skip dark particles */
  /*   fscanf(inputfile, "%d\n", &gid); */
  /* } */
  /* for(i=0;i<theader.nstar;i++){ // Add up Star Metals */
  /*   fscanf(inputfile, "%d\n", &gid); */
  /* } */
  /* fclose(inputfile); */
  /* printf("Reading .grp DONE.\n"); */
  return nskid_grps;
}

int read_grp(char *filename, int typecode)
{
  int i, nparts, dump;
  FILE *fp;
  if(!(fp = fopen(filename, "r")))
    {printf("Can't open .grp file! Force exit.\n"); exit(-1);}
  fscanf(fp, "%d\n", &nparts); // ONLY HERE: hid = Number of 
  if(typecode==0) nparts=theader.nsph;
  if(typecode==1) nparts=theader.nsph + theader.ndark;
  if(typecode==4) nparts=theader.nsph + theader.ndark + theader.nstar;
  galid = (int *)malloc(nparts*sizeof(int));
  i = 0;
  while(i<nparts){
    fscanf(fp, "%d\n", &dump);
    galid[i++] = dump;
    if(feof(fp)) break;
  }
  fclose(fp);
  for(i=0;i<nparts;i++)
    if(galid[i] < 0)
      galid[i] *= -1;
  fprintf(stdout, "---> Reading .grp done.\n");  
}

int read_sogrp(char *filename, int typecode)
{
  int i, nparts, dump;
  FILE *fp;
  if(!(fp = fopen(filename, "r")))
    {printf("Can't open .sogrp file! Force exit.\n"); exit(-1);}
  fscanf(fp, "%d\n", &nparts); // ONLY HERE: hid = Number of 
  if(typecode==0) nparts=theader.nsph;
  if(typecode==1) nparts=theader.nsph + theader.ndark;
  if(typecode==4) nparts=theader.nsph + theader.ndark + theader.nstar;
  sohid = (int *)malloc(nparts*sizeof(int));
  i = 0;
  while(i<nparts){
    fscanf(fp, "%d\n", &dump);
    sohid[i++] = dump;
    if(feof(fp)) break;
  }
  fclose(fp);
  for(i=0;i<nparts;i++)
    if(sohid[i] < 0)
      sohid[i] *= -1;
  fprintf(stdout, "---> Reading .sogrp done.\n");  
}

void get_snap_string(int snapnum, char *snapstr)
{
  char str[4];

  sprintf(str, "%d\0", snapnum);
  if (snapnum < 0)
    {printf("snapnum < 0! Force Exit!\n"); exit(-1);}
  else if (snapnum < 10)
    strcat(strcpy(snapstr,"00"), str);
  else if (snapnum < 100)
    strcat(strcpy(snapstr,"0"), str);
  else if (snapnum < 1000)
    strcat(strcpy(snapstr,""), str);
  else
    {printf("snapnum > 999! Force Exit!\n"); exit(-1);}
}

int headerinfo(struct tipsy_header *header)
{
  printf("*** Header information ***\n");
  printf("header.time: %f\n", header->time);
  printf("header.nbodies: %d\n", header->nbodies);
  printf("header.ndim: %d\n", header->ndim);
  printf("header.nsph: %d\n", header->nsph);
  printf("header.ndark: %d\n", header->ndark);
  printf("header.nstar: %d\n", header->nstar);
  printf("**************************\n");
}

void read_header(struct tipsy_header *head, FILE *ftipsy ) {
  fread((char *)&head->time, sizeof(head->time), 1, ftipsy);
  fread((char *)&head->nbodies, sizeof(head->nbodies), 1, ftipsy);
  fread((char *)&head->ndim, sizeof(head->ndim), 1, ftipsy);
  fread((char *)&head->nsph, sizeof(head->nsph), 1, ftipsy);
  fread((char *)&head->ndark, sizeof(head->ndark), 1, ftipsy);
  fread((char *)&head->nstar, sizeof(head->nstar), 1, ftipsy);
  //fread((char *)&head->pad, sizeof(head->pad), 1, ftipsy);
}
