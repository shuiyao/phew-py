#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "iontab.h"

#define XH 	0.76
#define KBOLTZ	1.381e-16
#define MHYDR	1.673e-24
#define CLIGHT	2.99792458e10
#define PI	3.14159265

double redshift;
char ionbkgd[15];
int nions;

typedef struct ionStruct {
  char name[10];
  float lambda,fraction,Xsec,atomwt,bsys,alpha;
  int Zcolumn;
} ionStruct;

extern ionStruct Ion[MAXIONS];
