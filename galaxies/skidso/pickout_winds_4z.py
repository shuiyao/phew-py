# Pick out all the winds that are ejected at z=0.1, 1.0, 2.0, 4.0
# Write new files winds.z0, winds.z1, winds.z2, winds.z4

import ioformat
import sys
import os

zs = [0.1, 1.0, 2.0, 4.0]
zstrs = ["z0", "z1", "z2", "z4"]

MODEL = "P50N288"

if(MODEL == "P25N144"):
    amin = ["0.909091", "0.5", "0.333333", "0.2"]
    amax = ["0.919091", "0.51", "0.343333", "0.21"]
    NCPU = 128
if(MODEL == "P50N288"):
    amin = ["0.909091", "0.5", "0.333333", "0.2"]
    amax = ["0.913091", "0.502", "0.335333", "0.202"]
    NCPU = 256
if(MODEL == "P50N576"):
    amin = ["0.909091", "0.5", "0.333333", "0.2"]
    amax = ["0.909291", "0.5002", "0.333533", "0.2002"]
    NCPU = 512

errormsg = "Usage: "+sys.argv[0]+" ncpu $WINDS"
if(len(sys.argv) != 3):
    print errormsg
    # sys.exit(1)
    FBASE = "/data002/shuiyao/data/p50n288fofIII/WINDS/"
    FBASE = "/data002/shuiyao/data/p25n144fofXIIIa/WINDS/"
else:
    NCPU = int(sys.argv[1])
    # First Wind have Tc([4]) = 1000 K
    # The tclock([2]) also provides information
    MODE = 0 # Prepare for the idlist
    FBASE = sys.argv[2]
    if(os.path.exists(FBASE) == 0):
        print errormsg
        print "Error: ", FBASE, "not found. Exit."
        sys.exit(1)
    print "Ncpu = ", NCPU

for zi in range(len(zs)):
    count = 0
    ofile = open(FBASE+"winds."+zstrs[zi], "w")
    print "Writting: ", FBASE+"winds."+zstrs[zi]
    for icpu in range(NCPU):
        ifile = open(FBASE+"winds."+str(icpu), "r")
        if(int(icpu/50) == icpu/50.):
            print "Doing File #", icpu
        for line in ifile:
            spt = line.split()
            if amin[zi] < spt[0] < amax[zi]:
                count += 1
                ofile.write(line)
    print "Total number of Winds: ", count
    ofile.close()
