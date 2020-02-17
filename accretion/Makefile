#adopted from tipsy-tools
CC = gcc
CFLAGS = -g -fPIC -O3
TOOLS = allstars
#OPT += -DMARK64
#OBJS = gadget3io.o
#LD_LIBRARY_PATH = /scratch/shuiyao/sci/gadget3io/
INCLS = tipsydefs.h /usr/include/stdio.h -I/scratch/shuiyao/sci/gadget3io/
LIBS = -L/scratch/shuiyao/sci/gadget3io/ -lm

#OPT = -DGALPARTS_INCLUDE_GAS # Find gas particles for galaxies?

all: $(TOOLS)

allstars: ${OBJS} allstars.o
	$(CC) $(CFLAGS) $(OPT) allstars.o $(OBJS) $(LIBS) -o allstars

allstars.o: allstars.c $(INCL)
	$(CC) $(CFLAGS) $(OPT) -c allstars.c $(LIBS)

clean:
	rm -f *.o
	rm -f $(TOOLS)


