/* Calculate the Mass fraction of a certain Ion */
/* fion ionid redshift n_h[cm^-3] T[K] */
/* IonID   0  1  2  3  4  5  6  7  8 */
/* IonName HI  HeII  CIII  CIV  OIV  OVI  NeVIII  MgII  SiIV */

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include "proto.h"
#include "allvars.h"
#include "iontab.h"
#include <unistd.h>

ionStruct Ion[MAXIONS];
float nhl[NHPTS],tl[TPTS];
float fractab[MAXIONS][NHPTS][TPTS];

int main(int argc,char **argv)
{
  /* double H_0; */
  /* float CosmicTime(); */
  double fracion;
  double nh, temp;
  int ionid;

  strcpy(ionbkgd,"HMQG");
#ifdef UVBKG_HM12
  strcpy(ionbkgd,"HM12");
#endif

  sscanf(argv[1],"%d",&ionid) ;  
  sscanf(argv[2],"%lg",&redshift) ;
  sscanf(argv[3],"%lg",&nh) ;
  sscanf(argv[4],"%lg",&temp) ;

  /* totMass = 0.30; //0.28; //0.238; //0.28; */
  /* Lambda = 0.70; //0.72; //0.762; //0.72; */
  /* Omega_b = 0.045; //0.046; //0.0418; //0.046; */
  /* H_0 = 70; //73; */
  /* h = 0.01*H_0; */

  //  cosmopar(CosmicTime(argv[1]));
#ifdef VERBOSE
  fprintf(stdout, "Initializing Ions ...\n");
#endif
  InitIons();
  // Should come before load_fraction_tables() to determine nions  
  
#ifdef VERBOSE
  fprintf(stdout, "Loading Ion Table ...\n");
#endif
  load_fraction_tables();
  
  fracion = IonFrac(temp, nh*MHYDR, ionid);
  /* temp in K */

  fprintf(stdout, "Z = %g, NH = %g, TEMPERATURE = %g\n", redshift, nh, temp);
  fprintf(stdout, "f(%s) = %g\n", Ion[ionid].name, fracion);

  return 0;
}

/* Read from file specions.dat:
   ion_name  rest_freq  osc_strength  atom_wt  solar_abundance  rel_abundance
 * solar abundance mass fractions are taken from Arnett, "Supernovae 
	and Nucleosynthesis" (1996) Table A.1 column 4,
	which in turn is from Anders & Grevesse 1989 and others.
 * rel_abundance is in dex vs.solar (not used in specgen)
*/
int InitIons()
{
  int i;
  char line[80],prefix[150],specionfilename[180];
  FILE *specfile;

  sprintf(prefix,"/home/shuiyao/temp/specexbin/ionfiles/");
  /* Default: specions_i9.dat */
  sprintf(specionfilename,"%sspecions_i9.dat",prefix);
  if( (specfile = fopen(specionfilename,"r")) == NULL ) specfile = fopen(specionfilename,"r");
  if( specfile == NULL ) {
    fprintf(stderr,"cannot find specion file anywhere\n");
    exit(-1);
  }

  i = 0;
  while( fgets(line,80,specfile) != NULL ) {
    if( strstr(line,"#") != NULL ) continue;
    if( i >= MAXIONS ) break;
    sscanf(line,"%10s %g %g %g %g %d %g",Ion[i].name,&Ion[i].lambda,&Ion[i].Xsec,&Ion[i].atomwt,&Ion[i].fraction,&Ion[i].Zcolumn,&Ion[i].alpha);
    i++;
  }
  nions = i;
  fclose(specfile);

#ifdef VERBOSE  
  fprintf(stderr,"Processing %d ions from specions.dat:\n",nions);
#endif  
  for( i=0; i<nions; i++ ) {
    Ion[i].bsys = sqrt(2.*KBOLTZ/(MHYDR*Ion[i].atomwt))/1.e5;
    Ion[i].Xsec *= 2.648e-2*Ion[i].lambda*1.e-13;
#ifdef VERBOSE      
    fprintf(stderr,"%5d %10s %12.6g %12.6g %10.5g %10.5g %10.5g % 3.1f % 2d\n",i,Ion[i].name,Ion[i].lambda,Ion[i].Xsec,Ion[i].atomwt,Ion[i].fraction,Ion[i].bsys*sqrt(1.e4),Ion[i].alpha,Ion[i].Zcolumn);
#endif    
  }
  return 0;
}
