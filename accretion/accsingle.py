# Track Accretion for a single simulation

import ioformat
import matplotlib.pyplot as plt
from scipy import sqrt, array, log10, linspace
from pylab import setp
from astroconst import pc, ac
from scipy import interpolate
import sys

TMAX_THRESH = 5.5
PLOT_MODE = 1 # Cold/Hot/Wind
# PLOT_MODE = 1 # Metal
idxlist = [1,3]
#idxlist = [0,2]

def tcosmic(a):
    if(a == 0): return 0
    else: return 13.7 * a / sqrt(0.30/(a**3) + 0.70)

class acc_halo:
    def __init__(self, mh):
        self.mh = mh
        self.cold = array([0.]*100)
        self.hot = array([0.]*100)
        self.wind = array([0.]*100)
        self.windret = array([0.]*100)
        self.zwind = array([0.]*100)
        # self.zprimo = array([0.]*100)
        self.zcold = array([0.]*100)
        self.zhot = array([0.]*100)
        self.dtime = array([0.]*100)        

ztab, ttab = ioformat.rcol("/home/shuiyao/code/Python/tcosmic.dat", [1,2], linestart=1)
tcosmic = interpolate.interp1d(ztab, ttab)

errormsg = "Usage: "+sys.argv[0]+" BoxSize[Mpc/h] $FBASE $simname $flag_phew"
if(len(sys.argv) != 5):
    print errormsg
    sys.exit(1)
else:
    BOXSIZE = float(sys.argv[1])
    FBASE = sys.argv[2]
    modelname = sys.argv[3]
    VOLUME = BOXSIZE ** 3
    fname = FBASE+modelname+"/SFRINFO/c.sfrinfo"
    foutname = "/scratch/shuiyao/scidata/newwind/acc."+modelname
    FFORMAT = int(sys.argv[4])
    print "modelname = ", modelname
    print "BoxSize = ", BOXSIZE, "Mpc/h"

zbins = linspace(0., 10., 100)
acchalos = [acc_halo(10.5)]
acchalos.append(acc_halo(11.0))
acchalos.append(acc_halo(11.5))
acchalos.append(acc_halo(12.0))

if(FFORMAT == 1):
    a, alast, tmax, mgal, mass, wmass, zmet = ioformat.rcol(fname, [0,2,3,4,5,6,7])
    for i in range(len(a)):
        mfof = mgal[i]
        if(mfof < 10.5): idx = 0
        elif(mfof < 11.0): idx = 1
        elif(mfof < 11.5): idx = 2
        else: idx = 3
        z = 1./a[i] - 1.
        zidx = int(z/0.1)
        if(zidx > 99): zidx = 99
        if(acchalos[idx].dtime[zidx] == 0.0):
            if(zidx==99): acchalos[idx].dtime[zidx] = tcosmic(zbins[zidx])
            else:
                t1, t2 = tcosmic(zbins[zidx+1]), tcosmic(zbins[zidx])
                acchalos[idx].dtime[zidx] = t2 - t1
        if(alast[i] == 0): # First Accretion
            if(tmax[i] < TMAX_THRESH):
                acchalos[idx].cold[zidx] += mass[i]
                acchalos[idx].zcold[zidx] += (mass[i] - wmass[i]) * zmet[i]
            else:
                acchalos[idx].hot[zidx] += mass[i]
                acchalos[idx].zhot[zidx] += (mass[i] - wmass[i]) * zmet[i]                      # acchalos[idx].zprimo[zidx] += mass[i] * zmet[i]
        elif(alast[i] < 0): # Wind re-accretion
            acchalos[idx].wind[zidx] += mass[i]
            acchalos[idx].zwind[zidx] += mass[i] * zmet[i]
    # ---------------- New Format ----------------
if(FFORMAT == 2):
    a, alast, tmax, mgal, mass, wmass, zmet = ioformat.rcol(fname, [0,2,3,4,5,6,7])
    # if fidx == 0: mass = array(mass) * ac.msolar * ac.msolar * 1.e20
    # if fidx == 0: wmass = array(wmass) * ac.msolar * ac.msolar * 1.e20
    for i in range(len(a)):
        mfof = mgal[i]
        if(mfof < 10.5): idx = 0
        elif(mfof < 11.0): idx = 1
        elif(mfof < 11.5): idx = 2
        else: idx = 3
        z = 1./a[i] - 1.
        zidx = int(z/0.1)
        if(zidx > 99): zidx = 99
        if(acchalos[idx].dtime[zidx] == 0.0):
            if(zidx==99): acchalos[idx].dtime[zidx] = tcosmic(zbins[zidx])
            else:
                t1, t2 = tcosmic(zbins[zidx+1]), tcosmic(zbins[zidx])
                acchalos[idx].dtime[zidx] = t2 - t1
        if(alast[i] == 0): # First Accretion
            if(wmass[i] > mass[i]): print "Wmass > Mass!", i
            if(tmax[i] < TMAX_THRESH):
                acchalos[idx].cold[zidx] += mass[i] - wmass[i]
                acchalos[idx].zcold[zidx] += (mass[i] - wmass[i]) * zmet[i]
            else:
                acchalos[idx].hot[zidx] += mass[i] - wmass[i]
                acchalos[idx].zhot[zidx] += (mass[i] - wmass[i]) * zmet[i]                     
            # acchalos[idx].zprimo[zidx] += (mass[i] - wmass[i]) * zmet[i] 
            acchalos[idx].wind[zidx] += wmass[i]
            # acchalos[idx].zwind[zidx] += wmass[i] * zmet[i]                
        elif(alast[i] < 0): # Wind re-accretion
            # if(wmass[i] / mass[i] < 0.9): print "Wmass != Mass!", i
            acchalos[idx].wind[zidx] += mass[i]
            acchalos[idx].zwind[zidx] += mass[i] * zmet[i]
            acchalos[idx].windret[zidx] += mass[i]
        else: # Stripped Re-accretion
            acchalos[idx].wind[zidx] += wmass[i]
            # acchalos[idx].zwind[zidx] += wmass[i] * zmet[i]

a3 = 1.0
for idx in range(4):
    acchalos[idx].cold /= (VOLUME * acchalos[idx].dtime * a3 * ac.msolar)
    acchalos[idx].hot /= (VOLUME * acchalos[idx].dtime * a3 * ac.msolar)
    acchalos[idx].wind /= (VOLUME * acchalos[idx].dtime * a3 * ac.msolar)
    if(FFORMAT == 2):
        acchalos[idx].windret /= (VOLUME * acchalos[idx].dtime * ac.msolar)
    acchalos[idx].zcold /= (VOLUME * acchalos[idx].dtime * ac.msolar)
    acchalos[idx].zhot /= (VOLUME * acchalos[idx].dtime * ac.msolar)
    acchalos[idx].zwind /= (VOLUME * acchalos[idx].dtime * ac.msolar)

suffix = [".105", ".110", ".115", ".inf"]    
for idx in range(4):
    fout = open(foutname+suffix[idx], "w")
    for i in range(len(zbins)):
        line = "% .3f % .5e % .5e % .5e % .5e % .5e % .5e % .5e\n" % \
               (zbins[i], acchalos[idx].cold[i], acchalos[idx].hot[i], \
                acchalos[idx].wind[i], acchalos[idx].windret[i], \
                acchalos[idx].zcold[i], acchalos[idx].zhot[i],
                acchalos[idx].zwind[i])
        fout.write(line)
    fout.close()
