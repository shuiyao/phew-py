import ioformat
import time
import sys
import os

# Select PhEW tracks from phews.? based on wid.dat
# use dictionary

# First Wind have Tc([4]) = 1000 K
# The tclock([2]) also provides information
errormsg = "Usage: "+sys.argv[0]+" ncpu redshift(0,1,2) $WINDS"
if(len(sys.argv) != 4):
    print errormsg
    sys.exit(1)
else:
    NCPU = int(sys.argv[1])
    REDSHIFT = float(sys.argv[2])
    FBASE = sys.argv[3]
    amaxs = ["1.0", "1.0", "1.0", "1.0"]
    if(REDSHIFT == 4.0):
        AMIN, AMAX, odir = "0.200", amaxs[0], FBASE+"z4/"
    elif(REDSHIFT == 2.0):
        AMIN, AMAX, odir = "0.333", amaxs[1], FBASE+"z2/"
    elif(REDSHIFT == 1.0):
        AMIN, AMAX, odir = "0.500", amaxs[2], FBASE+"z1/"
    elif(REDSHIFT == 0.2):
        AMIN, AMAX, odir = "0.833", amaxs[3], FBASE+"z0/"
    else:
        print errormsg
        sys.exit(1)
    if(os.path.exists(odir) == 0):
        print errormsg
        print "Error: ", FBASE, "not found. Exit."
        sys.exit(1)
    print "Ncpu = ", NCPU
    print "Directory: ", odir

widfile = odir+"wid.dat"
wids = ioformat.rcol(widfile, [0], [0])
print len(wids), "Wind IDs read..."

tstart = time.time()
tmid = tstart
for icpu in range(NCPU):
    print "Search File #", icpu, " | t = ", time.time() - tmid
    tmid = time.time()
    ifile = open(FBASE+"phews."+str(icpu), "r")
    ofile = open(odir+"phews."+str(icpu), "w")
    ifile.readline()
    for line in ifile:
        spt = line.split()
        if(float(spt[0]) > 1): continue # a > 1, never should happen
        wid = int(spt[1])
        if wid in wids:
            ofile.write(line)
    ifile.close()
    ofile.close()

print "Elapsed Time: ", time.time() - tstart, " sec"
