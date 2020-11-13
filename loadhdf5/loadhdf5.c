/* 180630 [CRITICAL] Bug Fix. Fix the mass and other properties for stars. */

#include <stdio.h>
#include "gadgetdefs.h"
#include "tipsydefs.h"
#include "loadhdf5.h"
#include <hdf5.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Unit Conversion for .tipsy */
/*   - Read in .hdf file, write it into .tipsy file. */
/*   - Usage: hdf52tipsy snap_* (no .hdf5) */

/* WARNING: */
/* Tipsy metallicity array contains only C, O, Si, Fe. Make sure the corresponding field is loaded from the HDF5 output. Default: [C, O, Si, Fe] = [2, 4, 7, 10] (NMETALS = 11) */

static long NumPart;

double unit_Time;
double unit_Length;
double unit_Density;
double unit_Mass;
double unit_Velocity;
double unit_Temp;

struct particle_data *P;
struct gadget_dump gheader;
struct gadget_dump h;

int flag_read_dark = 0;
int flag_read_star = 0;

#define MAX_NSKID_GRPS 50000

struct tipsy_header theader, soheader, galheader;
struct star_particle *gals, *sohalos;
struct skid_group *skid_grps;
struct fof_group *fofgrps;

int *sohid; /* Map from pidx -> HID */
int *galid; /* Map from pidx -> GID */
int *hosthid;

int Nfofgrps;

void cosmounits(void)
{
  double HubbleParam, BoxSize;
  HubbleParam = gheader.HubbleParam;
  BoxSize = gheader.BoxSize;
  unit_Time=sqrt(8*M_PI/3)*CM_PER_MPC/(100*HubbleParam*1.e5);
  unit_Density=1.8791E-29*HubbleParam*HubbleParam;
  unit_Length=BoxSize*CM_PER_MPC*1.e-3;
  unit_Mass=unit_Density*unit_Length*unit_Length*unit_Length/(HubbleParam*HubbleParam);
  unit_Velocity=unit_Length/unit_Time;
  unit_Temp = pow(UNIT_V, 2);
  fprintf(stdout, "unit_Length = %g\n", unit_Length);
  fprintf(stdout, "unit_Mass = %g\n", unit_Mass);
  fprintf(stdout, "unit_Velocity = %g\n", unit_Velocity);
  // = pow(UNIT_L, 2) / pow((UNIT_L/Unit_V), 2); 
  return;
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

int load_hdf5(char *snap){

  char infile[256];
  FILE *fp;
  int multipart = 0;

  int allocate_memory();
  
  sprintf(infile,"%s.hdf5",snap);
  if(!(fp=fopen(infile,"r"))){
    fprintf(stdout, "%s not found.\n", infile);
    sprintf(infile,"%s.0.hdf5",snap);
    multipart = 1;
    if(!(fp=fopen(infile,"r"))){
      fprintf(stderr,"Error opening file '%s' \n",infile);
      exit(0);
    }
  }
  fclose(fp);

  long i, cnt;
  int j, k;
  long noffset;
  hid_t hdf5_file, hdf5_headergrp, hdf5_attribute, hdf5_grp, hdf5_dataset;
  /* int ngas  = 1; // Why? */
  long ngas  = 0; 
  long ndark  = 0; 
  long nstar  = 0; 
  float *posvel, *single, *metals;
  int *intsingle;

  hdf5_file = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT);
  hdf5_headergrp = H5Gopen1(hdf5_file, "/Header");
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"NumFilesPerSnapshot");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &h.num_files);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Redshift");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &h.redshift);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"BoxSize");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &h.BoxSize);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Omega0");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &h.Omega0);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"OmegaLambda");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &h.OmegaLambda);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Time");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &h.time);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"HubbleParam");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &h.HubbleParam);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"NumPart_Total");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, h.npartTotal);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"NumPart_ThisFile");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, h.npart);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"MassTable");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, h.mass);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Flag_Metals");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &h.flag_metals);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Flag_StellarAge");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &h.flag_stellarage);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Flag_Sfr");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &h.flag_sfr);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Flag_Cooling");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &h.flag_cooling);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Flag_Feedback");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &h.flag_feedback);
  H5Aclose(hdf5_attribute);

  /* Missing several fields including flag_sfr, flag_feedback, nparttotal, etc. */

  H5Gclose(hdf5_headergrp);
  H5Fclose(hdf5_file);

  NumPart = h.npartTotal[0] + h.npartTotal[1] + h.npartTotal[4];

  h.flag_metals = NMETALSHDF5;

  gheader = h;

  fprintf(stderr,"NumFiles = %d\n", h.num_files);  
  fprintf(stderr,"NumPart = %ld\n", NumPart);
  fprintf(stderr,"Time = %g; Redshift = %g\n", h.time, h.redshift);  
  
  allocate_memory();

  for(k=0; k<h.num_files; k++){
    if(multipart)
      sprintf(infile,"%s.%d.hdf5",snap,k);

    hdf5_file = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT);
    if(multipart){
      hdf5_headergrp = H5Gopen1(hdf5_file, "/Header");
      hdf5_attribute = H5Aopen_name(hdf5_headergrp,"NumPart_ThisFile");
      H5Aread(hdf5_attribute, H5T_NATIVE_INT, h.npart);
      H5Aclose(hdf5_attribute);
      hdf5_attribute = H5Aopen_name(hdf5_headergrp,"MassTable");
      H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, h.mass);
      H5Aclose(hdf5_attribute);
      H5Gclose(hdf5_headergrp);
    }

    //GAS
    fprintf(stdout, "Reading GAS for file #%d\n", k);    
    if(!(posvel = (float *)malloc(sizeof(float)*h.npart[0]*3))){
      fprintf(stderr, "Failed to allocate memory for posvel\n");
      exit(-1);
    }
    if(!(single = (float *)malloc(sizeof(float)*h.npart[0]))){
      fprintf(stderr, "Failed to allocate memory for single\n");
      exit(-1);
    }
    if(!(metals = (float *)malloc(sizeof(float)*h.npart[0]*h.flag_metals))){
      fprintf(stderr, "Failed to allocate memory for metals\n");
      exit(-1);
    }
#ifdef PHEW_TRACK_INFO
    float *trackinfo;
    if(!(trackinfo = (float *)malloc(sizeof(float)*h.npart[0]*4))){
      fprintf(stderr, "Failed to allocate memory for trackinfo\n");
      exit(-1);
    }
#endif    
    if(!(intsingle = (int *)malloc(sizeof(int)*h.npart[0]))){
      fprintf(stderr, "Failed to allocate memory for intsingle\n");
      exit(-1);
    }
    
    hdf5_grp = H5Gopen1(hdf5_file, "/PartType0");
    
    hdf5_dataset = H5Dopen1(hdf5_grp, "Coordinates");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      for(j=0; j<3; j++)
	P[i].Pos[j] = posvel[cnt*3 + j];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "Velocities");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      for(j=0; j<3; j++)
	P[i].Vel[j] = posvel[cnt*3 + j];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "ParticleIDs");
    H5Dread(hdf5_dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, intsingle);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].ID = intsingle[cnt];
      cnt += 1;
    }

    if(h.mass[0] == 0 && h.npart[0] > 0){
      hdf5_dataset = H5Dopen1(hdf5_grp, "Masses");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
      H5Dclose(hdf5_dataset);
      for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
	P[i].Mass = single[cnt];
	cnt += 1;
      }
    }
    else{
      for(i=ngas; i<h.npart[0]+ngas; i++)
	P[i].Mass = h.mass[0];
    }
    
    hdf5_dataset = H5Dopen1(hdf5_grp, "InternalEnergy");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Temp = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "Density");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Rho = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "ElectronAbundance");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Ne = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "NeutralHydrogenAbundance");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Nh = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "SmoothingLength");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Hsml = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "StarFormationRate");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++) {
      P[i].Sfr = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "TemperatureMax");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Tmax = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "DelayTime");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].DelayTime = single[cnt];
      cnt += 1;
    }

    if(H5Lexists(hdf5_grp,"FractionH2",H5P_DEFAULT)==0){
      for(i=ngas; i<h.npart[0]+ngas; i++){
	P[i].fH2 = 0.0;
      }
    }
    else{
      hdf5_dataset = H5Dopen1(hdf5_grp, "FractionH2");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
      H5Dclose(hdf5_dataset);
      for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
	P[i].fH2 = single[cnt];
	cnt += 1;
      }
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "Metallicity");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, metals);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      for(j=0; j<h.flag_metals; j++)
	P[i].metal[j] = metals[cnt*h.flag_metals + j];
      cnt += 1;
    }

#ifdef PHEW_EXTRA_OUTPUT    
    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWKey");
    H5Dread(hdf5_dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, intsingle);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Key = intsingle[cnt];
      cnt += 1;
    }
    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWMcloud");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Mcloud = single[cnt];
      cnt += 1;
    }
    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWWindMass");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].WindMass = single[cnt];
      cnt += 1;
    }
    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWLastSFTime");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].LastSFTime = single[cnt];
      cnt += 1;
    }
    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWVinit");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      P[i].Vinit = single[cnt];
      cnt += 1;
    }
#endif
#ifdef PHEW_TRACK_INFO
    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWTrackInfo");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, trackinfo);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<h.npart[0]+ngas; i++){
      for(j=0; j<4; j++)
	P[i].TrackInfo[j] = trackinfo[cnt*4+j];
      cnt += 1;
    }
#endif    

    H5Gclose(hdf5_grp);
#ifdef PHEW_TRACK_INFO
    free(trackinfo);
#endif    
    free(metals);
    free(single);
    free(posvel);
    ngas += h.npart[0];

    // DARK
    if(flag_read_dark){
      fprintf(stdout, "Reading DARK from file #%d\n", k);    
      noffset = gheader.npartTotal[0];
      if(!(posvel = (float *)malloc(sizeof(float)*h.npart[1]*3))){
	fprintf(stderr, "Failed to allocate memory for posvel\n");
	exit(-1);
      }
      if(!(single = (float *)malloc(sizeof(float)*h.npart[1]))){
	fprintf(stderr, "Failed to allocate memory for single\n");
	exit(-1);
      }
      if(!(intsingle = (int *)malloc(sizeof(int)*h.npart[1]))){
	fprintf(stderr, "Failed to allocate memory for intsingle\n");
	exit(-1);
      }
    
      hdf5_grp = H5Gopen1(hdf5_file, "/PartType1");

      hdf5_dataset = H5Dopen1(hdf5_grp, "Coordinates");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel);
      H5Dclose(hdf5_dataset);
      for(i=noffset+ndark, cnt=0; i<noffset+h.npart[1]+ndark; i++){
	for(j=0; j<3; j++)
	  P[i].Pos[j] = posvel[cnt*3 + j];
	cnt += 1;
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "Velocities");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel);
      H5Dclose(hdf5_dataset);
      for(i=noffset+ndark, cnt=0; i<noffset+h.npart[1]+ndark; i++){
	for(j=0; j<3; j++)
	  P[i].Vel[j] = posvel[cnt*3 + j];
	cnt += 1;
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "ParticleIDs");
      H5Dread(hdf5_dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, intsingle);
      H5Dclose(hdf5_dataset);
      for(i=noffset+ndark, cnt=0; i<noffset+h.npart[1]+ndark; i++){
	P[i].ID = intsingle[cnt];
	cnt += 1;
      }

      if(h.mass[1] == 0 && h.npart[1] > 0){
	hdf5_dataset = H5Dopen1(hdf5_grp, "Masses");
	H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
	H5Dclose(hdf5_dataset);
	for(i=noffset+ndark, cnt=0; i<noffset+h.npart[1]+ndark; i++){
	  P[i].Mass = single[cnt];
	  cnt += 1;
	}
      }
      else
	for(i=noffset+ndark; i<noffset+h.npart[1]+ndark; i++)
	  P[i].Mass = h.mass[1];

      H5Gclose(hdf5_grp);
      free(single);
      free(posvel);
      ndark += h.npart[1];
    } // flag_read_dark

    if(flag_read_star){
      fprintf(stdout, "Reading STAR from file #%d\n", k);
      // STAR
      noffset = gheader.npartTotal[0] + gheader.npartTotal[1];
      if(!(posvel = (float *)malloc(sizeof(float)*h.npart[4]*3))){
	fprintf(stderr, "Failed to allocate memory for posvel\n");
	exit(-1);
      }
      if(!(single = (float *)malloc(sizeof(float)*h.npart[4]))){
	fprintf(stderr, "Failed to allocate memory for single\n");
	exit(-1);
      }
      if(!(metals = (float *)malloc(sizeof(float)*h.npart[4]*h.flag_metals))){
	fprintf(stderr, "Failed to allocate memory for single\n");
	exit(-1);
      }
      if(!(intsingle = (int *)malloc(sizeof(int)*h.npart[4]))){
	fprintf(stderr, "Failed to allocate memory for intsingle\n");
	exit(-1);
      }
    
      hdf5_grp = H5Gopen1(hdf5_file, "/PartType4");

      hdf5_dataset = H5Dopen1(hdf5_grp, "Coordinates");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel);
      H5Dclose(hdf5_dataset);
      for(i=noffset+nstar, cnt=0; i<noffset+h.npart[4]+nstar; i++){
	for(j=0; j<3; j++)
	  P[i].Pos[j] = posvel[cnt*3 + j];
	cnt += 1;
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "Velocities");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel);
      H5Dclose(hdf5_dataset);
      for(i=noffset+nstar, cnt=0; i<noffset+h.npart[4]+nstar; i++){
	for(j=0; j<3; j++)
	  P[i].Vel[j] = posvel[cnt*3 + j];
	cnt += 1;
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "ParticleIDs");
      H5Dread(hdf5_dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, intsingle);
      H5Dclose(hdf5_dataset);
      for(i=noffset+nstar, cnt=0; i<noffset+h.npart[4]+nstar; i++){
	P[i].ID = intsingle[cnt];
	cnt += 1;
      }

      if(h.mass[4] == 0 && h.npart[4] > 0){
	hdf5_dataset = H5Dopen1(hdf5_grp, "Masses");
	H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
	H5Dclose(hdf5_dataset);
	for(i=noffset+nstar, cnt=0; i<noffset+h.npart[4]+nstar; i++){
	  P[i].Mass = single[cnt];
	  cnt += 1;
	}
      }
      else
	for(i=noffset+nstar; i<noffset+h.npart[4]+nstar; i++)
	  P[i].Mass = h.mass[4];

      /* For Star particles */
      hdf5_dataset = H5Dopen1(hdf5_grp, "StellarFormationTime");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
      H5Dclose(hdf5_dataset);
      for(i=noffset+nstar, cnt=0; i<noffset+h.npart[4]+nstar; i++){
	P[i].Sfr = single[cnt];
	cnt += 1;
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "TemperatureMax");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
      H5Dclose(hdf5_dataset);
      for(i=noffset+nstar, cnt=0; i<noffset+h.npart[4]+nstar; i++){
	P[i].Tmax = single[cnt];
	cnt += 1;
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "Metallicity");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, metals);
      H5Dclose(hdf5_dataset);
      for(i=noffset+nstar, cnt=0; i<noffset+h.npart[4]+nstar; i++){
	for(j=0; j<h.flag_metals; j++)
	  P[i].metal[j] = metals[cnt*h.flag_metals + j];
	cnt += 1;
      }
    
      H5Gclose(hdf5_grp);
      free(posvel);
      free(single);
      free(metals);
      free(intsingle);
      nstar += h.npart[4];
    } // flag_read_star

    H5Fclose(hdf5_file);
    fprintf(stderr, "File: %d ngas = %d(%5.3f)\n", k, h.npart[0],
	    (float)(h.npart[0])/(float)(NumPart));
    fprintf(stderr, "File: %d ndark = %d(%5.3f)\n", k, h.npart[1],
	    (float)(h.npart[1])/(float)(NumPart));
    fprintf(stderr, "File: %d nstar = %d(%5.3f)\n", k, h.npart[4],
	    (float)(h.npart[4])/(float)(NumPart));
  }
  return 0;
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

void read_sopar(char *filename)
{
  int i;
  FILE *fp;
  int gid, hid;
  
  hosthid = (int *) malloc(MAX_NSKID_GRPS * sizeof(int));
  i = 0;
  if(!(fp = fopen(filename, "r"))){
    fprintf(stderr, "Can't open .par file! Force exit.\n"); exit(-1);}
  while(1){
    fscanf(fp, "%d %d", &gid, &hid);
    hosthid[gid] = hid;
    if(feof(fp)) break;
    i++;
  }
  fclose(fp);
}

int read_skid(int snapnum)
{
  FILE *inputfile;
  char skidname[200], skidbasename[200], grpname[200];
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
  /* for(i=0;i<nparts;i++) */
  /*   if(sohid[i] < 0) */
  /*     sohid[i] *= -1; */
  fprintf(stdout, "---> Reading .sogrp done.\n");  
}

int allocate_memory(void)
{
  fprintf(stdout, "Allocating %6.3f GB Memory for %ld particles.\n",
	  NumPart * sizeof(struct particle_data) / (1024. * 1024. * 1024.),
	  NumPart);

  if(!(P=malloc(NumPart * sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }

  fprintf(stdout, "Memory allocated.\n");
  /* P--;    */
  /* start with offset 1 */
  return 0;
}
