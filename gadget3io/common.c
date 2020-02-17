#include <stdio.h>
#include <string.h>
#include "tipsydefs.h"
#include "proto.h"

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

int headerinfo(struct dump *header)
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

void read_header(struct dump *head, FILE *ftipsy ) {
  fread((char *)&head->time, sizeof(head->time), 1, ftipsy);
  fread((char *)&head->nbodies, sizeof(head->nbodies), 1, ftipsy);
  fread((char *)&head->ndim, sizeof(head->ndim), 1, ftipsy);
  fread((char *)&head->nsph, sizeof(head->nsph), 1, ftipsy);
  fread((char *)&head->ndark, sizeof(head->ndark), 1, ftipsy);
  fread((char *)&head->nstar, sizeof(head->nstar), 1, ftipsy);
  //fread((char *)&head->pad, sizeof(head->pad), 1, ftipsy);
}
