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
    # acosmic(tcosmic(amin) + 1.5e9) Track for 1.5 Gyr
    amaxs = ["0.318264", "0.434506", "0.593577", "0.934432"]
    if(REDSHIFT == 4.0):
        AMIN, AMAX, odir = "0.200", amaxs[0], FBASE+"z4/"
    elif(REDSHIFT == 2.0):
        AMIN, AMAX, odir = "0.333333", amaxs[1], FBASE+"z2/"
    elif(REDSHIFT == 1.0):
        AMIN, AMAX, odir = "0.500", amaxs[2], FBASE+"z1/"
    elif(REDSHIFT == 0.2):
        AMIN, AMAX, odir = "0.833333", amaxs[3], FBASE+"z0/"
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
wids = set(wids) # Enabling fast search with wid in wids
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
        if(spt[0] > "1"): continue # a > 1, never should happen
        if(spt[0] < AMIN): continue # a < amin, Not selected
        if(spt[0] > AMAX): continue # a > amax (1.5 Gyr after amin)        
        wid = int(spt[1])
        if wid in wids: # <----- The time consuming line.
            ofile.write(line)
    ifile.close()
    ofile.close()

print "Elapsed Time: ", time.time() - tstart, " sec"
