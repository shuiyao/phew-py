import ioformat
import time
import sys
import os

# First Wind have Tc([4]) = 1000 K
# The tclock([2]) also provides information
errormsg = "Usage: "+sys.argv[0]+" ncpu redshift(0,1,2) $WINDS"
if(len(sys.argv) != 4):
    print errormsg
    sys.exit(1)
else:
    NCPU = int(sys.argv[1])
    # First Wind have Tc([4]) = 1000 K
    # The tclock([2]) also provides information
    MODE = 0 # Prepare for the idlist
    REDSHIFT = float(sys.argv[2])
    FBASE = sys.argv[3]
    if(REDSHIFT == 4.0):
        AMIN, AMAX, odir = "0.200", "0.210", FBASE+"z4/"
    elif(REDSHIFT == 2.0):
        AMIN, AMAX, odir = "0.333", "0.433", FBASE+"z2/"
    elif(REDSHIFT == 1.0):
        AMIN, AMAX, odir = "0.500", "0.600", FBASE+"z1/"
    elif(REDSHIFT == 0.2):
        AMIN, AMAX, odir = "0.833", "1.000", FBASE+"z0/"
    else:
        print errormsg
        sys.exit(1)
    if(os.path.exists(odir) == 0):
        print errormsg
        print "Error: ", FBASE, "not found. Exit."
        sys.exit(1)
    print "Ncpu = ", NCPU
    print "Redshift = ", REDSHIFT
    print "amin, amax = ", AMIN, AMAX
    print "Directory: ", odir

MODE = 0 # Prepare for the idlist
widfile = odir+"wid.dat"

wids = ioformat.rcol(widfile, [0], [0])
idmax = max(wids)
idmask = [0]*(idmax+1)
for wid in wids: idmask[wid] = 1
print len(wids), "Wind IDs read..."

tstart = time.time()
tmid = tstart
for icpu in range(NCPU):
    print "Search File #", icpu, " | t = ", time.time() - tmid
    tmid = time.time()
    # ifile = open("analytic_track."+str(icpu), "r")
    # ofile = open(odir+"analytic_track."+str(icpu), "w")
    ifile = open(FBASE+"track."+str(icpu), "r")
    ofile = open(odir+"track."+str(icpu), "w")    
    for line in ifile:
        spt = line.split()
        if(float(spt[0]) > 1): continue
        wid = int(spt[1])
        if wid <= idmax:
            if idmask[wid] == 1:
                if(AMIN < spt[0] < AMAX):
                    ofile.write(line)
    ifile.close()
    ofile.close()

print "Elapsed Time: ", time.time() - tstart, " sec"
