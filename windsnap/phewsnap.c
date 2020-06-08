#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include "tipsydefs.h"
#include "proto.h"
#ifdef IONS
#include "iontab.h"
#endif

#define NUMOFNODES_X 288
#define NUMOFNODES_Y 288
#define UNIT_M 433697.735404 // 47963.49995379916 [P12]
#define UNIT_L 25000.0
#define HUBBLEPARAM 0.7
#define XH 0.76
#define MHYDR 1.673e-24
#define MCLOUD0 1.e5 // Initial mass of each cloudlet

/* #define LSIZE (600.0/12000.0) */
/* #define UNIT_M 47963.49995379916 */
/* #define UNIT_L 12000.0 */

/* Call: phewsnap fbin faux modelname faw zmin zmax */

/* Output: $SCIDATA/modelname/, halo_$num.grid; halo_$num.winds */
/* halo_$num.grid: */
/* Mvir, Rvir, xcen, ycen, zcen */
/* xcells, ycells */
/* (Integrated) Mass, Mass(O), (MassHI, MassOVI) */

typedef struct NodesStruct {
  float x;
  float y;
  float m;
  float mo;
#ifdef IONS
  float mhi;
  float movi;
#endif  
} NodesStruct;
NodesStruct *nodes_array;

int ncells_x;
int ncells_y;
char infilename[200], outfilename[200];
double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Xcen, Ycen, Zcen;
double Vxcen, Vycen, Vzcen;
double Mvir, Rvir;
int Ntot = 0, Nwind = 0;
#ifdef IONS
double a3inv, unit_Density;
double fHI, fOVI;
ionStruct Ion[MAXIONS];
#endif

struct dump header, soheader;
struct gas_particle *gp;
struct dark_particle *dp;
struct star_particle *sp;
struct aux_gas_data *aux_gp;
struct aux_star_data *aux_sp;
struct star_particle *sohalo;
struct aw_gas_data *aw_gp;

int LoadGrid(int nx, int ny)
{
  int i, nnodes, ix, iy;
  float dx, dy, x, y;
  int ncells;
  ncells = nx * ny;
  nodes_array = (struct NodesStruct *)
    malloc(ncells*sizeof(*nodes_array));
  dx = (Xmax - Xmin) / nx;
  dy = (Ymax - Ymin) / ny;
  for(ix=0;ix<nx;ix++) {
    x = Xmin + dx * ix;
    for (iy=0;iy<ny;iy++){
      i = ix * ny + iy;
      nodes_array[i].x = x;
      nodes_array[i].y = Ymin + dy * iy;
      nodes_array[i].m = 0.;
      nodes_array[i].mo = 0.;
    }
  }
  printf("Grid Loaded.\n");
  return 1;
}

int GridCount(float x, float y, float z, float m, int flag, float meto)
{
  int ix, iy, i;
  if(z > Zmin && z < Zmax){
    ix = (int)((x - Xmin) / (Xmax - Xmin) * ncells_x);
    if (ix < ncells_x && ix > 0){
      iy = (int)((y - Ymin) / (Ymax - Ymin) * ncells_y);
      if (iy < ncells_y && iy > 0){
	if(flag == 0){
	i = ix * ncells_y + iy;
	nodes_array[i].m += m;
	nodes_array[i].mo += m * meto;
#ifdef IONS
	nodes_array[i].mhi += m * XH * fHI;
	nodes_array[i].movi += m * meto * fOVI;
#endif	
	Ntot ++;
	return 0;
	}
	else{ // A wind
	  Nwind ++;
	  Ntot ++;
	  return 1;
	}
      }
    }
  }
  return 0;
}

void WriteGrid(char *outgridfilename, char *outwindfilename, int flag)
{
  FILE *Fout;
  int ix, iy, i;
  double vdotr;
  Fout = fopen(outgridfilename, "w");
  fprintf(Fout, "%g %g %g %g\n", UNIT_L, Xcen, Ycen, Zcen);  
  fprintf(Fout, "%d %d\n", ncells_x, ncells_y);
  /* for(ix=0;ix<ncells_x;ix++) {fprintf(Fout, "%g\n", nodes_array[ix*ncells_y].x);} */
  /* for(iy=0;iy<ncells_y;iy++) {fprintf(Fout, "%g\n", nodes_array[iy].y);} */
  for(ix=0;ix<ncells_x;ix++) {
    for (iy=0;iy<ncells_y;iy++){
      i = ix * ncells_y + iy;
      fprintf(Fout, "%g %g", 
	      nodes_array[i].m, nodes_array[i].mo);
	      /* nodes_array[i].met * 1.28571); */
#ifdef IONS
      fprintf(Fout, " %g %g\n",
	      nodes_array[i].mhi, nodes_array[i].movi);
#else
      fprintf(Fout, "\n");       
#endif      
    }
  }
  fclose(Fout);

  Fout = fopen(outwindfilename, "w");
  for(i = 0; i < header.nsph; i++){
    if(aw_gp[i].wind_flag == 1){ /* It's marked as wind. */
      // In Gadget3, P.Mass*SphP.n.NumNgb = 4./3.*PI*SphP.Hsml**3*SphP.d.Density
      fprintf(Fout, "%g %g %g %g %g %g %g %g %g %g %g %d\n",
	      i,
	      gp[i].pos[0] - Xcen,
	      gp[i].pos[1] - Ycen,
	      gp[i].pos[2] - Zcen,
	      gp[i].mass * UNIT_M * 1.e10,
	      aux_gp[i].metal[1] * gp[i].mass * UNIT_M * 1.e10,
	      gp[i].hsmooth * UNIT_L,
	      aw_gp[i].rcloud * UNIT_L,
	      aw_gp[i].temp,
	      gp[i].mass * UNIT_M * 1.e10 / aw_gp[i].mcloud / MCLOUD0, // Clumping factor
	      aw_gp[i].mcloud, aw_gp[i].lastsftime,
	      aw_gp[i].wind_flag);
    }
  }
  fclose(Fout);
}

int read_tipsy_binary(FILE *fp, int typecode)
{
#ifdef MARK64
  printf("I am reading for 64 machine!\n");
  read_header(&header, fp);
#else
  printf("I am reading for 32 machine!\n");
  if(fread((char *)&header,sizeof(header),1,fp) == 0)
    {printf("BREAK!\n"); exit(-1);}
#endif
  //  headerinfo(&header);
  printf("*** Header Information ***\n");
  printf("header.time: %f\n", header.time);
  printf("header.nbodies: %d\n", header.nbodies);
  printf("header.ndim: %d\n", header.ndim);
  printf("header.nsph: %d\n", header.nsph);
  printf("header.ndark: %d\n", header.ndark);
  printf("header.nstar: %d\n", header.nstar);
  printf("**************************\n");
  if(header.nsph != 0) {
    gp = (struct gas_particle *)
      malloc(header.nsph*sizeof(*gp));
    if(gp == NULL) {
      printf("<sorry, no memory for gas particles, master>\n") ;
      return -1;
    }
  }
  if(header.ndark != 0 && typecode >= 1) {
    dp = (struct dark_particle *)
      malloc(header.ndark*sizeof(*dp));
    if(dp == NULL) {
      printf("<sorry, no memory for dark particles, master>\n") ;
      return -1;
    }
  }
  if(header.nstar != 0 && typecode == 4) {
    sp = (struct star_particle *)
      malloc(header.nstar*sizeof(*sp));
    if(sp == NULL) {
      printf("<sorry, no memory for star particles, master>\n") ;
      return -1;
    }
  }
  fread((char *)gp,sizeof(struct gas_particle),	header.nsph,fp) ;
  if(typecode >= 1)
    fread((char *)dp,sizeof(struct dark_particle), header.ndark,fp) ;
  if(typecode == 4)
    fread((char *)sp,sizeof(struct star_particle), header.nstar,fp) ;
  fclose(fp);
  printf("---> Reading .bin done.\n");
  fprintf(stdout,"read time %f\n",header.time);
  return 1;
}

int read_tipsy_aux(FILE *fp, int typecode)
{
/* Read in auxiliary particle data file */
  aux_gp = (struct aux_gas_data *) malloc(header.nsph*sizeof(*aux_gp));
  fread(aux_gp,sizeof(struct aux_gas_data),header.nsph,fp);
  if(typecode == 4){
    aux_sp = (struct aux_star_data *) malloc(header.nstar*sizeof(*aux_sp));
    fread(aux_sp,sizeof(*aux_sp),header.nstar,fp);
  }
  fclose(fp);
  fprintf(stdout, "---> Reading .aux done.\n");  
  return 1;
}

int read_tipsy_aw(FILE *fp)
{
/* Read in auxiliary particle data file */
  aw_gp = (struct aw_gas_data *) malloc(header.nsph*sizeof(*aw_gp));
  fread(aw_gp,sizeof(struct aw_gas_data),header.nsph,fp);
  fclose(fp);
  fprintf(stdout, "---> Reading .aw done.\n");  
  return 1;
}

int read_sogtp(FILE *fp)
{
#ifdef MARK64
  printf("I am reading for 64 machine!\n");
  read_header(&soheader, fp);
#else
  printf("I am reading for 32 machine!\n");
  if(fread((char *)&soheader,sizeof(struct dump)-4,1,fp) == 0) // PAD HERE
    {printf("EMPTY SOGTP FILE, BREAK!\n"); exit(-1);}
#endif
/*   headerinfo(&header); */
  printf("*** Sogtp Information ***\n");
  printf("soheader.nstar: %d\n", soheader.nstar);
  printf("**************************\n");

  if(soheader.nsph != 0) {fprintf(stdout, "BAD: NSPH != 0\n"); exit(-1);}
  if(soheader.ndark != 0) {fprintf(stdout, "BAD: NDARK != 0\n"); exit(-1);}
  if(soheader.nstar != 0) {
    sohalo = (struct star_particle *)
      malloc(soheader.nstar*sizeof(*sohalo));
    if(sohalo == NULL) {
      printf("<sorry, no memory for star particles, master>\n") ;
      return -1;
    }
  }
  fread((char *)sohalo,sizeof(struct star_particle),
  	soheader.nstar,fp) ;
  fclose(fp);
  printf("---> Reading .sogtp done.\n");
  return 1;
}

int main(int argc, char **argv)
{
  /* Calling Sequence: */
  /* phewsnap binfile, auxfile, modelname, awfile */
  
  char outgridfilename[500], outwindfilename[500];
  char modelname[500];
  char dummy[500];
  FILE *Fbin, *Faux, *Faw;
  int i, halonum;
  int awflag = 0;

  strcpy(modelname, argv[3]);  

  if( (Fbin = fopen(argv[1],"r")) == NULL) {
    fprintf(stderr,"cannot open %s\n",argv[1]);
    exit(-1);
  }
  if( (Faux = fopen(argv[2],"r")) == NULL) {
    fprintf(stderr,"cannot open %s\n",argv[2]);
    exit(-1);
  }
  if(argc >= 5){
    if( (Faw = fopen(argv[4],"r")) == NULL) {
      fprintf(stderr,"cannot open %s\n",argv[4]);
      exit(-1);
    }
    awflag = 1;
  }
  if(argc >= 6){
    Zmin = atof(argv[5]);
    Zmax = atof(argv[6]);
  } else {
    Zmin = -0.5;
    Zmax = 0.5;
  }

  read_tipsy_binary(Fbin, 0);/* OUTPUT: struct (*_particle) *gp, *dp, *sp */
  read_tipsy_aux(Faux, 0);/* OUTPUT: struct aux_gas_data *aux_gp */
  if(argc >= 5) read_tipsy_aw(Faw);

#ifdef IONS
  double rho;
  unit_Density = 1.87e-29 * HUBBLEPARAM * HUBBLEPARAM;
  a3inv = 1. / (header.time * header.time * header.time);
  InitIons(1./header.time - 1.);  // Including load_fraction_tables()
  /* fHI = IonFrac(1.e3, 10.0*MHYDR, 0); // test */
#endif  

  strcat(strcpy(outgridfilename, "/scratch/shuiyao/scidata/windsnap/"), modelname);
  strcat(strcpy(outwindfilename, outgridfilename), "/halo_");
  strcat(outwindfilename, argv[3]);
  strcat(outwindfilename, ".wind");    
  strcat(outgridfilename, "/halo_");
  strcat(outgridfilename, argv[3]);
  strcat(outgridfilename, ".grid");

  fprintf(stdout, "Output files as:\n%s\n%s\n",
	  outgridfilename, outwindfilename);

  Xcen = 0.0;
  Ycen = 0.0;
  Xmin = -0.5;
  Xmax = 0.5;
  Ymin = -0.5;
  Ymax = 0.5;
  /* Zmin = -0.5; */
  /* Zmax = 0.5; */
  Zcen = (Zmin + Zmax) / 2.0;  

  fprintf(stdout, "Box Range: [%5.3f, %5.3f], [%5.3f, %5.3f], [%5.3f, %5.3f]\n",
	  Xmin, Xmax, Ymin, Ymax, Zmin, Zmax);

  ncells_x = NUMOFNODES_X;
  ncells_y = NUMOFNODES_Y;

  LoadGrid(ncells_x, ncells_y);

  for(i=0; i<header.nsph; i++){
#ifdef IONS
    fHI = 0.0;
    fOVI = 0.0;
    if(aux_gp[i].delaytime <= 0){
      rho = (gp[i].rho * HUBBLEPARAM * HUBBLEPARAM * XH) * unit_Density * a3inv;
      fHI = IonFrac(gp[i].temp, rho, 0);
      fOVI = IonFrac(gp[i].temp, rho, 5);
      if(gp[i].temp < 1.e5) fOVI = 0.0; // Only collisionally ionized OVI
    }
#endif    
    if(GridCount(gp[i].pos[0], gp[i].pos[1], gp[i].pos[2],
		 gp[i].mass, aw_gp[i].wind_flag, aux_gp[i].metal[1])) /* C, O, Si, Fe; */
      {} // Marked as wind
  }
  
  fprintf(stdout, "Grid Calculated: Nwind/Ntot = %d/%d\n", Nwind, Ntot);

  WriteGrid(outgridfilename, outwindfilename, awflag);
  free(nodes_array);
}
