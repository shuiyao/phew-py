OPT     +=  -DDO9IONS
#OPT	+=  -DVERBOSE
#OPT	+=  -DNONEQUTIL

CC= gcc
FC= f77
CLINK=gcc
FLINK=f77
CFLAGS= ${OPT} -O2 -g -Wall 
FFLAGS= -O
CLIB= -lm #-lhdf5
FLIB= 

all: fion

OBJS= fion.o ionfrac.o

fion:	$(OBJS) Makefile
	$(CLINK) $(CFLAGS) -o fion $(OBJS) $(CLIB)

clean:
	rm -f *.o
	rm -f fion

