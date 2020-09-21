import matplotlib.pyplot as plt
from astroconst import pc, ac
from numpy import log10, genfromtxt, array, inf, linspace
import h5py
import ioformat
# from scipy.interpolate import interp1d
from numpy import polyfit, polyval
import config_mpl
import bin1d

# Based on readhdf5.py

BOXSIZE = 50000
UNIT_MASS = 3469581.88
HPARAM = 0.7
XH = 0.76
XHE = (1.0 - XH) / (4.0 * XH)
UNIT_T = 1.e10 * (2./3.) * 1.6726e-24 / 1.3806e-16

# UNIT_MASS = 3469581.88 * (12./50.) ** 3

#modelname = "l50n288-phewoff"
#modelname = "l50n288-phew-m5"
modelname = "p50n288fiducial"
fbase = "/nas/astro-th/shuiyao/"
snapname = fbase+modelname+"/snapshot_108.hdf5"
galname = fbase+modelname+"/gal_z108.stat"
soname = fbase+modelname+"/so_z108.sovcirc"
gasname = fbase+modelname+"/so_z108.mgas"

boundsvalue=0.90
def plotmedian(x0, y0, ax, nbins=20, clr="blue", alphavalue=0.4, verbose=False):
    x, y = bin1d.subsample(x0, y0, nonzero=True)
    s = bin1d.bin1d(x, y, nbins=nbins, bounds=boundsvalue)
    xline, yline = s.cen, s.median
    uline, lline = s.ubound, s.lbound
    if(verbose == True):
        print ("xline and yline: ")
        print (len(x))
        print (xline, yline)
    for i in range(len(yline))[1:-1]:
        if(yline[i] == 0):
            yline[i] = 0.5 * (yline[i-1] + yline[i+1])
            uline[i] = 0.5 * (uline[i-1] + uline[i+1])
            lline[i] = 0.5 * (lline[i-1] + lline[i+1])            
    p, = ax.plot(xline, yline, "-", color=clr) # Draw line for all models    
    pf, = plt.fill(xline+xline[::-1], uline+lline[::-1], color=clr, alpha=alphavalue)

def polyfit_mstar_to_msub():
    halos = genfromtxt(soname, names=True)
    gals = genfromtxt(galname)
    mstars = []
    for gal in gals: mstars.append(log10(gal[4] * UNIT_MASS * 1.e10 / HPARAM))
    mstars.append(-inf)
    mstars = array(mstars)
    hostmask = (halos['Msub'] > 0) & (mstars > 0)
    mstars = mstars[hostmask]
    msubs = log10(halos[hostmask]['Msub'] / HPARAM)
    # mfits = interp1d(msubs[::10], mstars[::10], kind="quadratic")
    mask = (msubs > 11.0) & (mstars > 8.0)
    msubs = msubs[mask]
    mstars = mstars[mask]
    coefs = polyfit(mstars, msubs, 2)
    # x = linspace(8.0, 11.0, 100)
    # y = polyval(coefs, x)
    return coefs

coefs = polyfit_mstar_to_msub()
x = array([9.0, 9.5, 10.0, 10.5, 11.0])
y = polyval(coefs, x)
halos = genfromtxt(gasname, names=True)
halos = halos[halos['Msub'] > 0]
fcold = 10.0 ** (halos['Mcold'] - halos['Msub'])
fhot = 10.0 ** (halos['Mhot'] - halos['Msub'])
fism = 10.0 ** (halos['Mism'] - halos['Msub'])
step = 2

fig, ax1 = plt.subplots(1, 1, figsize=(8,7))
# ax1.plot(halos['Msub'][::step], fcold[::step], "b.", markersize=4)
# ax1.plot(halos['Msub'][::step], fhot[::step], "r.", markersize=4)
# ax1.plot(halos['Msub'][::step], fism[::step], "g.", markersize=4)
plotmedian(halos['Msub'], fcold, ax1, nbins=15, clr="blue", alphavalue=0.4, verbose=False)
plotmedian(halos['Msub'], fhot, ax1, nbins=15, clr="red", alphavalue=0.4, verbose=False)
plotmedian(halos['Msub'], fism, ax1, nbins=15, clr="orange", alphavalue=0.4, verbose=False)

ax1.set_xlim(11.0, 13.0)
ax1.set_xticks([11.0, 11.5, 12.0, 12.5, 13.0])
ax1.set_ylim(0.0, 0.2)
ax1.set_xlabel(r"$M_{vir}$")
ax1.plot([11.0, 13.0], [0.15, 0.15], "k--")
ax2 = ax1.twiny()
ax2.set_xticks(y.tolist())
ax2.set_xticklabels(["9.0", "9.5", "10.0", "10.5", "11.0"])
ax2.set_xlabel(r"$M_*$")
ax2.set_xlim(11.0, 13.0)
ax2.figure.canvas.draw()

from pltastro import legend
lgd = legend.legend()
# lgd.addLine(("Cold", 'blue', '-', 1))
# lgd.addLine(("Hot", 'red', '-', 1))
# lgd.addLine(("ISM", 'orange', '-', 1))
lgd.addPatch(("Cold", 'blue', 'black'))
lgd.addPatch(("Hot", 'red', 'black'))
lgd.addPatch(("ISM", 'orange', 'black'))
lgd.loc = "upper right"
lgd.fontsize=12
lgd.draw()
plt.savefig("/home/shuiyao_umass_edu/figures/tmp.pdf")

plt.show()

# plt.plot(msubs[::3], mstars[::3], "b.", alpha=0.4)
# plt.plot(x, y, "r-")


print ("DONE")
