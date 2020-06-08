#adopted from tipsy-tools
CC = gcc
CFLAGS = -O2 -g
TOOLS = windsnap phewsnap
OBJSPHEW = phewsnap.o common.o ionfrac.o
OBJSWIND = windsnap.o common.o ionfrac.o
OPT += -DIONS
#OPT += -DMARK64
#OPT += -DEXCLUDE_WIND_PARTICLES # Whether we do not add winds to column density (windsnap ONLY). By default phewsnap does NOT use wind particles because they are too small 

all: $(TOOLS) 

windsnap: $(OBJSWIND)
	$(CC) $(CFLAGS) $(OPT) -o windsnap $(OBJSWIND) $(LIBS) -lm

windsnap.o: windsnap.c /usr/include/stdio.h tipsydefs.h proto.h
	$(CC) $(CFLAGS) $(OPT) -c windsnap.c $(LIBS) -lm

phewsnap: $(OBJSPHEW)
	$(CC) $(CFLAGS) $(OPT) -o phewsnap $(OBJSPHEW) $(LIBS) -lm

phewsnap.o: phewsnap.c /usr/include/stdio.h tipsydefs.h proto.h
	$(CC) $(CFLAGS) $(OPT) -c phewsnap.c $(LIBS) -lm

ionfrac.o: ionfrac.c /usr/include/stdio.h proto.h iontab.h
	$(CC) $(CFLAGS) $(OPT) -c ionfrac.c $(LIBS) -lm

common.o: common.c /usr/include/stdio.h tipsydefs.h proto.h
	$(CC) $(CFLAGS) $(OPT) -c common.c $(LIBS)

clean:
	rm -f *.o
	rm -f $(TOOLS)


