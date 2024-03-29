import ioformat
import sys
import os


errormsg = "Usage: "+sys.argv[0]+" ncpu redshift(0,1,2) $WINDS"
if(len(sys.argv) != 4):
    print (errormsg)
    sys.exit(1)
else:
    NCPU = int(sys.argv[1])
    # First Wind have Tc([4]) = 1000 K
    # The tclock([2]) also provides information
    MODE = 0 # Prepare for the idlist
    REDSHIFT = float(sys.argv[2])
    amaxs = ["0.22", "0.35", "0.52", "0.85"]
    FBASE = sys.argv[3]
    if(REDSHIFT == 4.0):
        AMIN, AMAX, odir, zstr = "0.200", amaxs[0], FBASE, "z4"
    elif(REDSHIFT == 2.0):
        AMIN, AMAX, odir, zstr = "0.333333", amaxs[1], FBASE, "z2"
    elif(REDSHIFT == 1.0):
        AMIN, AMAX, odir, zstr = "0.500", amaxs[2], FBASE, "z1"
    elif(REDSHIFT == 0.2):
        AMIN, AMAX, odir, zstr = "0.833333", amaxs[3], FBASE, "z0"
    else:
        print (errormsg)
        sys.exit(1)
    if(os.path.exists(FBASE) == 0):
        print (errormsg)
        print ("Error: ", FBASE, "not found. Exit.")
        sys.exit(1)
    if(os.path.exists(odir) == 0):    
        os.mkdir(odir)
    print ("Ncpu = ", NCPU)
    print ("Redshift = ", REDSHIFT)
    print ("amin, amax = ", AMIN, AMAX)
    print ("Directory: ", odir)

ids = []
ofile = open(odir+"winds."+zstr, "w")    
for icpu in range(NCPU):
    ifile = open(FBASE+"winds."+str(icpu), "r")
    print ("Doing File #", icpu)
    for line in ifile:
        spt = line.split()
        if AMIN < spt[0] < AMAX:
            # The Ids are PhEWKey 
            ids.append(spt[1])
            ofile.write(line)
    ifile.close()
ofile.close()
