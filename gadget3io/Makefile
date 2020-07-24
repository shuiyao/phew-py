#adopted from tipsy-tools
CC = gcc
CFLAGS = -g
TOOLS = allstars
#TOOLS = allstars get_halo_particles
#OPT += -DMARK64
OBJS = gadget3io.o
INCL = gadget3io.h tipsydefs.h /usr/include/stdio.h
LIBS = -lm

#OPT = -DGALPARTS_INCLUDE_GAS # Find gas particles for galaxies?

all: $(TOOLS)

gadget3io: gadget3io.o
	$(CC) $(CFLAGS) $(OPT) -o gadget3io gadget3io.o $(LIBS)

gadget3io.o: gadget3io.c $(INCL)
	$(CC) $(CFLAGS) $(OPT) -c gadget3io.c $(LIBS)

cpairing: cpairing.o ${OBJS} 
	$(CC) $(CFLAGS) $(OPT) cpairing.o $(OBJS) $(LIBS) -o cpairing

cpairing.o: cpairing.c $(INCL)
	$(CC) $(CFLAGS) $(OPT) -c cpairing.c $(LIBS)

galparts: ${OBJS} galparts.o
	$(CC) $(CFLAGS) $(OPT) galparts.o $(OBJS) $(LIBS) -o galparts

galparts.o: galparts.c $(INCL)
	$(CC) $(CFLAGS) $(OPT) -c galparts.c $(LIBS)

allstars: ${OBJS} allstars.o
	$(CC) $(CFLAGS) $(OPT) allstars.o $(OBJS) $(LIBS) -o allstars

allstars.o: allstars.c $(INCL)
	$(CC) $(CFLAGS) $(OPT) -c allstars.c $(LIBS)

get_halo_particles: ${OBJS} get_halo_particles.o
	$(CC) $(CFLAGS) $(OPT) get_halo_particles.o $(OBJS) $(LIBS) -o get_halo_particles

get_halo_particles.o: get_halo_particles.c $(INCL)
	$(CC) $(CFLAGS) $(OPT) -c get_halo_particles.c $(LIBS)

# skidsopid: ${OBJS} skidsopid.o
# 	$(CC) $(CFLAGS) $(OPT) skidsopid.o $(OBJS) $(LIBS) -o skidsopid

# skidsopid.o: skidsopid.c $(INCL)
# 	$(CC) $(CFLAGS) $(OPT) -c skidsopid.c $(LIBS)

clean:
	rm -f *.o
	rm -f $(TOOLS)


