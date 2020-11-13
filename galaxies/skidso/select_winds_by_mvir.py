import ioformat
import sys
import os

errormsg = "Usage: "+sys.argv[0]+" ncpu redshift(0,1,2) $WINDS"
if(len(sys.argv) != 4):
    print errormsg
    MODE = 0
    NCPU = 128
    FBASE = "/scratch/shuiyao/data/p25n144rwlIII/WINDS/"
    AMIN, AMAX, odir, num = "0.500", "0.510", FBASE+"z1/", "005"
    fpidgrps = FBASE+"../pid_z"+num+".grp"
    # sys.exit(1)
else:
    NCPU = int(sys.argv[1])
    # First Wind have Tc([4]) = 1000 K
    # The tclock([2]) also provides information
    MODE = 0 # Prepare for the idlist
    REDSHIFT = float(sys.argv[2])
    FBASE = sys.argv[3]
    if(REDSHIFT == 2.0):
        AMIN, AMAX, odir, num = "0.333", "0.335", FBASE+"z2/", "003"
    elif(REDSHIFT == 1.0):
        AMIN, AMAX, odir, num = "0.500", "0.502", FBASE+"z1/", "005"
    elif(REDSHIFT == 0.2):
        AMIN, AMAX, odir, num = "0.833", "0.835", FBASE+"z0/", "100"
    else:
        print errormsg
        sys.exit(1)
    if(os.path.exists(FBASE) == 0):
        print errormsg
        print "Error: ", FBASE, "not found. Exit."
        sys.exit(1)
    if(os.path.exists(odir) == 0):    
        os.mkdir(odir)
    print "Find Winds According to Mvir ---> "
    print "Ncpu = ", NCPU
    print "Redshift = ", REDSHIFT
    print "amin, amax = ", AMIN, AMAX
    print "Directory: ", odir
    fpidgrps = FBASE+"../pid_z"+num+".grp"

mvirs = ioformat.rcol(fpidgrps, [3])
# pidx, hid, Mgal, Mvir, Msub
ids = [[], [], [], [], []]
for icpu in range(NCPU):
    ifile = open(FBASE+"winds."+str(icpu), "r")
    print "Doing File #", icpu
    for line in ifile:
        spt = line.split()
        if AMIN < spt[0] < AMAX:
            pid = int(spt[1])
            mfof = float(spt[4])
            mvir = mvirs[pid-1]
            if(abs(mvir - mfof) < 1.0):
                if(10.3 < mvir < 10.7): ids[0].append(spt[1])
                elif(10.9 < mvir < 11.1): ids[1].append(spt[1])
                elif(11.4 < mvir < 11.6): ids[2].append(spt[1])
                elif(11.9 < mvir < 12.1): ids[3].append(spt[1])
                elif(12.5 < mvir < 13.0): ids[4].append(spt[1])                                

for k in range(5): print len(ids[k])," Wind Record Found!"
k = 0
for mh in ["105", "110", "115", "120", "125"]:
    ofile = open(odir+"wid"+mh+".dat", "w")
    for wid in ids[k]:
        ofile.write(wid+"\n")
    ofile.close()
    k = k + 1
