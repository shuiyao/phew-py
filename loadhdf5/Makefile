#HDF5INCL = -I/Users/bob/local/hdf5/include
# HDF5INCL = -I/usr/include/hdf5/serial
# HDF5LIB  = -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a -lpthread -lz -ldl

# HDF5INCL = -I/home/shuiyao/hdf5/include
HDF5LIB = -lhdf5 -lz

OPT += -DH5_USE_16_API
#OPT += -DNMETALS=11 # How many species are tracked in the HDF5 "Metallicity"  block? If not set correctly, code returns segmentation fault after the output "Reading GAS for file #%d"

OPT += -DPHEW_EXTRA_OUTPUT
OPT += -DPHEW_TRACK_INFO

OPT += -DNOSPHSKIP

CC= gcc
CLINK=gcc
#CFLAGS= -O2 -g -Wall $(OPT) $(HDF5INCL)
CFLAGS= -g -fPIC $(OPT) $(HDF5INCL)
CLIB= -lm $(HDF5LIB)

INCL = loadhdf5.h gadgetdefs.h tipsydefs.h

all: rhot get_particles get_halo_particles

OBJS= loadhdf5.o

loadhdf5.o: loadhdf5.c $(INCL)
	$(CC) $(CFLAGS) -c loadhdf5.c $(CLIB)

get_particles: get_particles.o $(OBJS) 
	 $(CLINK) $(CFLAGS) get_particles.o $(OBJS) $(CLIB) -o get_particles

get_particles.o: get_particles.c $(OBJS) $(INCL)
	 $(CLINK) $(CFLAGS) -c get_particles.c $(CLIB)

rhot: rhot.o $(OBJS) 
	 $(CLINK) $(CFLAGS) rhot.o $(OBJS) $(CLIB) -o rhot

rhot.o: rhot.c $(OBJS) $(INCL)
	 $(CLINK) $(CFLAGS) -c rhot.c $(CLIB)

get_halo_particles: $(OBJS) get_halo_particles.o 
	 $(CLINK) $(CFLAGS) get_halo_particles.o $(OBJS) $(CLIB) -o get_halo_particles

get_halo_particles.o: get_halo_particles.c $(INCL)
	 $(CLINK) $(CFLAGS) -c get_halo_particles.c $(CLIB)

clean:
	rm -f ./*.o
	rm -f loadhdf5

