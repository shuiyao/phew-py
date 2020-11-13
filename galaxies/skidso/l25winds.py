#import matplotlib.pyplot as plt
import ioformat
from numpy import logspace, array, linspace, sqrt, histogram, log10, histogram2d, pi, median
from astroconst import pc, ac
import sys
import os

DEBUG = False

# VinitVc(wps, 4, foutname="output.dat")
# write_wind_features(wps)
errormsg = "Usage: "+sys.argv[0]+" redshift(0,1,2,4) $fbase $simname"
if(len(sys.argv) != 4):
    print errormsg
    sys.exit(1)
else:
    MODE = 0 # Prepare for the idlist
    REDSHIFT = float(sys.argv[1])
    FBASE = sys.argv[2]
    modelname = sys.argv[3]
    ascale = 1./(REDSHIFT+1.)
    if(REDSHIFT == 4.0):
        AMIN, AMAX, snapnum, subwfolder = 0.20, 0.22, "033", "z4"        
    elif(REDSHIFT == 2.0):
        AMIN, AMAX, snapnum, subwfolder = 0.333333, 0.353333, "058", "z2"
    elif(REDSHIFT == 1.0):
        AMIN, AMAX, snapnum, subwfolder = 0.50, 0.52, "078", "z1"
    elif(REDSHIFT == 0.2):
        AMIN, AMAX, snapnum, subwfolder = 0.82, 0.846667, "108", "z0"
    else:
        print errormsg
        sys.exit(1)
    if(os.path.exists(FBASE) == 0):
        print errormsg
        print "Error: ", FBASE, "not found. Exit."
        sys.exit(1)
    print "Redshift = ", REDSHIFT
    print "amin, amax = ", AMIN, AMAX
    print "Directory: ", FBASE+modelname

database = FBASE+modelname
grpbase = FBASE+modelname
outputbase = "/data002/shuiyao/scripts/windsinfo/"
fname = database+"/WINDS/"+subwfolder+"/sorted.tracks"
sogrp = grpbase+"/so_z"+snapnum+".sogrp"
sovcirc = grpbase+"/so_z"+snapnum+".sovcirc"
stat = grpbase+"/gal_z"+snapnum+".stat"
fid = database+"/snapshot_"+snapnum+".id"
fwind = outputbase+modelname+"/windsinfo."+subwfolder
#unit_m = 5995.4321936
unit_m = 433697.735404
Mgasp = 9.3e7

if (MODE == 0):
    mvir, rvir = ioformat.rcol(sovcirc, [1,2], linestart=1)
    mstar = ioformat.rcol(stat, [4])
    pidl = ioformat.rcol(fid, [0], [0])
    hid = ioformat.rcol(sogrp, [0], [0], linestart=1)
    pidx = [-1] * (len(hid)+1)
    for i in range(len(pidl)):
        if(pidx[pidl[i]] == -1): # Ensure that only gas PID is valid
            pidx[pidl[i]] = i # pidx maps PID to idx
    mvir = array(mvir) / 0.7
    rvir = array(rvir) / 0.7 * ascale
    mstar = array(mstar) / 0.7 * unit_m * 1.e10
    # vc = sqrt(pc.G * mvir * ac.msolar / (rvir * ac.kpc)) / 1.e5
    # plt.plot(mvir, vc, "b.", markersize=2)
    # plt.xscale("log")
    # plt.show()

#wpbins = select()
print "compiled."

OMEGA0 = 0.30
OMEGA_LAMBDA = 0.70
def tcosmic(a):
    z = 1./(a) - 1.
    return 13.7 * 1.e9 * a / sqrt(OMEGA0/a**3 + OMEGA_LAMBDA)
def overdensity(z):
    x = OMEGA_LAMBDA / (OMEGA0 * (1.+z)**3)
    x = 18. * pi * pi * (1. + 0.4093*x**0.9052)
    return x

class windp():
    def __init__(self, pid):
        self.a = []
        self.pid = pid
        self.mvir = 0.
        self.rvir = 0.
        self.mstar = 0.
        self.tclock = []
        self.dr = []
        self.dr2 = []
        self.mh = []
        self.dv = []
        self.dt = []        

def select():
    a, pid, tclock, mh, dr, dv, rho, T = \
    ioformat.rcol(fname, [0,1,2,3,4,5,6,7], [1])

    idmax = 0
    i, j = 0, 0
    wps = []
    for i in range(len(mh)):
        if(pid[i] > idmax): # Find the vrel bin, The first appearance of this ID
            ainit = a[i]
            idmax = pid[i]
            if(AMIN <= a[i] <= AMAX):
                wps.append(windp(pid[i]))
                wps[-1].a.append(a[i])
                t0 = tcosmic(a[i])
                idx = pidx[pid[i]]
                if(abs(hid[idx]) != 0):
                    wps[-1].mvir = mvir[abs(hid[idx])-1]
                    wps[-1].rvir = rvir[abs(hid[idx])-1]
                    wps[-1].mstar = mstar[abs(hid[idx])-1]
        if(AMIN <= ainit <= AMAX):
            wps[-1].dt.append((tcosmic(a[i])-t0)/1.e6) # Myr
            wps[-1].dr.append(dr[i])
            wps[-1].dv.append(dv[i])
            wps[-1].mh.append(mh[i])
    return wps

def write_wind_features(wps, foutname=fwind):
    fout = open(foutname, "w")
    fout.write("#ID Mvir Rvir Vinit V25 Rret\n")
    count, writecount = 0, 0
    for w in wps:
        if(max(w.dr) < 1000. and 0.1 < w.mvir*0.7/10.**w.mh[0] < 20. and w.rvir > 0):
            mratio = w.mvir*0.7/10.**w.mh[0]
            vc, v25, vinit = 0., 0., 0.            
            x, y = [w.dr[0]], [w.dv[0]]
            ddr, ddt = [], []
            i25, r25 = 0, -1
            vel25 = -1
            ireturn = 0
            count += 1
            vc = sqrt(pc.G * w.mvir * ac.msolar / (w.rvir * ac.kpc)) / 1.e5
            for i in range(len(w.dr))[1:]:
                x.append(abs(w.dr[i]))
                y.append(w.dv[i])
                ddt.append(w.dt[i]-w.dt[i-1])
                ddr.append(x[i]-x[i-1])
                if(i25 == 0 and abs(w.dr[i]) > w.rvir * 0.25 and w.dt[i]<400.): # WARNING: dt < 400.
                    i25 = i
                if(ireturn == 0 and ((ddr[i-1] < 0. and w.dv[i] < vc) or ddt[i-1]>30.)):
                    ireturn = i
                    rreturn = x[i]
            if(ireturn == 0):
                ireturn = -1
                rreturn = -1
            if(max(ddr) < 100. and min(ddr) > -100.):
                writecount += 1
                if(i25 > 0):
                    v25 = y[i25]
                else:
                    v25 = -1
                vinit = y[0]
                # fout.write("#ID, Mvir, Vc, Vinit, V25\n")            
                outstr = str(w.pid)+" "+str(log10(w.mvir))+" "+str(w.rvir)
                outstr += " "+str(vinit)+" "+str(v25)+" "+str(rreturn)
                outstr += "\n"
                fout.write(outstr)
    print writecount, count, len(wps)
    fout.close()

if(DEBUG == False):
    print "Now select wind particles..."
    wps = select()
    print "Now write: ", fwind
    write_wind_features(wps)
    
