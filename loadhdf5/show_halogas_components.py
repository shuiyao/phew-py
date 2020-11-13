import matplotlib.pyplot as plt
from astroconst import pc, ac
from numpy import log10, genfromtxt, array, inf, linspace
import h5py
import ioformat
# from scipy.interpolate import interp1d
from numpy import polyfit, polyval
import config_mpl
import bin1d
from matplotlib.pyplot import setp

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
models = ["p50n288fiducial", "l50n288-phewoff", "l50n288-phew-m5"]
modellgd = ["Gadget3", "GIZMO-PhEWOff", "GIZMO-PhEW"]
lstyles = [":", "--", "-"]
# fbase = "/nas/astro-th-nas/shuiyao/"
fbase = "/nas/astro-th/shuiyao/"
# snapstr = "108"
# snapname = fbase+modelname+"/snapshot_"+snapstr+".hdf5"
# galname = fbase+modelname+"/gal_z"+snapstr+".stat"
# soname = fbase+modelname+"/so_z"+snapstr+".sovcirc"
# gasname = fbase+modelname+"/so_z"+snapstr+".mbar.T5"

boundsvalue=0.68
def plotmedian(x0, y0, ax, nbins=20, clr="blue", lstyle="-", alphavalue=0.4, fill=True, verbose=False):
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
    p, = ax.plot(xline, yline, linestyle=lstyle, color=clr) # Draw line for all models
    if(fill == True):
        pf, = ax.fill(xline+xline[::-1], uline+lline[::-1], color=clr, alpha=alphavalue)

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

fbase = "/nas/astro-th-nas/shuiyao/"
#fbase = "/nas/astro-th/shuiyao/"

fig, axs = plt.subplots(1, 2, figsize=(9,6))
ax1 = axs[0]
snapstr = "108"
suffix = ".T80k"

for mi, modelname in enumerate(models):
    snapname = fbase+modelname+"/snapshot_"+snapstr+".hdf5"
    galname = fbase+modelname+"/gal_z"+snapstr+".stat"
    soname = fbase+modelname+"/so_z"+snapstr+".sovcirc"
    gasname = fbase+modelname+"/so_z"+snapstr+".mbar"+suffix

    coefs = polyfit_mstar_to_msub()
    x = array([9.0, 9.5, 10.0, 10.5, 11.0])
    y = polyval(coefs, x)
    halos = genfromtxt(gasname, names=True)
    halos = halos[halos['Msub'] > 0]
    fcold = 10.0 ** (halos['Mcold'] - halos['Msub'])
    fhot = 10.0 ** (halos['Mhot'] - halos['Msub'])
    fism = 10.0 ** (halos['Mism'] - halos['Msub'])
    ftot = fcold + fhot + fism + 10.0 ** (halos['Mstar'] - halos['Msub'])
    fcold = fcold / 0.15
    fhot = fhot / 0.15
    fism = fism / 0.15
    ftot = ftot / 0.15    

    step = 2

    # ax1.plot(halos['Msub'][::step], fcold[::step], "b.", markersize=4)
    # ax1.plot(halos['Msub'][::step], fhot[::step], "r.", markersize=4)
    # ax1.plot(halos['Msub'][::step], fism[::step], "g.", markersize=4)
    if(modelname == "l50n288-phew-m5"):
        plotmedian(halos['Msub'], fcold, ax1, nbins=15, lstyle=lstyles[mi], clr="blue", alphavalue=0.4, verbose=False)
        plotmedian(halos['Msub'], fhot, ax1, nbins=15, lstyle=lstyles[mi], clr="red", alphavalue=0.4, verbose=False)
        plotmedian(halos['Msub'], fism, ax1, nbins=15, lstyle=lstyles[mi], clr="orange", alphavalue=0.4, verbose=False)
        plotmedian(halos['Msub'], ftot, ax1, nbins=15, lstyle=lstyles[mi], clr="grey", alphavalue=0.4, verbose=False)
    else:
        plotmedian(halos['Msub'], fcold, ax1, nbins=15, lstyle=lstyles[mi], clr="blue", fill=False, verbose=False)
        plotmedian(halos['Msub'], fhot, ax1, nbins=15, lstyle=lstyles[mi], clr="red", fill=False, verbose=False)
        plotmedian(halos['Msub'], fism, ax1, nbins=15, lstyle=lstyles[mi], clr="orange", fill=False, verbose=False)
        plotmedian(halos['Msub'], ftot, ax1, nbins=15, lstyle=lstyles[mi], clr="grey", fill=False, verbose=False)

ax1.set_xlim(11.0, 13.0)
ax1.set_xticks([11.0, 11.5, 12.0, 12.5])
# ax1.set_ylim(0.0, 0.16)
# ax1.set_yticks([0.0, 0.04, 0.08, 0.12, 0.16])
ax1.set_ylim(0.0, 1.0)
# ax1.set_yticks([0.0, 0.04, 0.08, 0.12, 0.16])
ax1.set_xlabel(r"$\log(M_{vir})$")
# ax1.plot([11.0, 13.0], [0.15, 0.15], "k--")
ax2 = ax1.twiny()
ax2.set_xticks(y.tolist())
ax2.set_xticklabels(["9.0", "9.5", "10.0", "10.5", "11.0"])
ax2.set_xlabel(r"$\log(M_*)$")
ax2.set_xlim(11.0, 13.0)
ax2.figure.canvas.draw()
ax1.set_ylabel(r"$M_{gas}/M_{vir}$")
# ax1.set_title("z = 0")
ax1.set_title("log(T/K) = 4.9")

from pltastro import legend
lgd = legend.legend(ax1)
# lgd.addLine(("Cold", 'blue', '-', 1))
# lgd.addLine(("Hot", 'red', '-', 1))
# lgd.addLine(("ISM", 'orange', '-', 1))
lgd.addPatch(("Cold", 'blue', 'black'))
lgd.addPatch(("Hot", 'red', 'black'))
lgd.addPatch(("ISM", 'orange', 'black'))
lgd.addPatch(("All", 'grey', 'black'))
lgd.loc = "upper right"
lgd.fontsize=12
lgd.draw()

ax1b = axs[1]
snapstr = "108"
#snapstr = "108"
suffix = ".T6"
for mi, modelname in enumerate(models):
    snapname = fbase+modelname+"/snapshot_"+snapstr+".hdf5"
    galname = fbase+modelname+"/gal_z"+snapstr+".stat"
    soname = fbase+modelname+"/so_z"+snapstr+".sovcirc"
    gasname = fbase+modelname+"/so_z"+snapstr+".mbar"+suffix

    coefs = polyfit_mstar_to_msub()
    x = array([9.0, 9.5, 10.0, 10.5, 11.0])
    y = polyval(coefs, x)
    halos = genfromtxt(gasname, names=True)
    halos = halos[halos['Msub'] > 0]
    fcold = 10.0 ** (halos['Mcold'] - halos['Msub'])
    fhot = 10.0 ** (halos['Mhot'] - halos['Msub'])
    fism = 10.0 ** (halos['Mism'] - halos['Msub'])
    ftot = fcold + fhot + fism + 10.0 ** (halos['Mstar'] - halos['Msub'])
    fcold = fcold / 0.15
    fhot = fhot / 0.15
    fism = fism / 0.15
    ftot = ftot / 0.15    
    if(modelname == "l50n288-phew-m5"):
        plotmedian(halos['Msub'], fcold, ax1b, nbins=15, lstyle=lstyles[mi], clr="blue", alphavalue=0.4, verbose=False)
        plotmedian(halos['Msub'], fhot, ax1b, nbins=15, lstyle=lstyles[mi], clr="red", alphavalue=0.4, verbose=False)
        plotmedian(halos['Msub'], fism, ax1b, nbins=15, lstyle=lstyles[mi], clr="orange", alphavalue=0.4, verbose=False)
        plotmedian(halos['Msub'], ftot, ax1b, nbins=15, lstyle=lstyles[mi], clr="grey", alphavalue=0.4, verbose=False)
    else:
        plotmedian(halos['Msub'], fcold, ax1b, nbins=15, lstyle=lstyles[mi], clr="blue", fill=False, verbose=False)
        plotmedian(halos['Msub'], fhot, ax1b, nbins=15, lstyle=lstyles[mi], clr="red", fill=False, verbose=False)
        plotmedian(halos['Msub'], fism, ax1b, nbins=15, lstyle=lstyles[mi], clr="orange", fill=False, verbose=False)
        plotmedian(halos['Msub'], ftot, ax1b, nbins=15, lstyle=lstyles[mi], clr="grey", fill=False, verbose=False)

ax1b.set_xlim(11.0, 13.0)
ax1b.set_xticks([11.0, 11.5, 12.0, 12.5])
# ax1b.set_ylim(0.0, 0.16)
# ax1b.set_yticks([0.0, 0.04, 0.08, 0.12, 0.16])
ax1b.set_ylim(0.0, 1.0)
# ax1b.set_ylabel(r"$M_{gas}/M_{vir}$")
setp(ax1b.get_yticklabels(), visible=False)

ax1b.set_xlabel(r"$\log(M_{vir})$")
# ax1b.plot([11.0, 13.0], [0.15, 0.15], "k--")
ax2b = ax1b.twiny()
ax2b.set_xticks(y.tolist())
ax2b.set_xticklabels(["9.0", "9.5", "10.0", "10.5", "11.0"])
ax2b.set_xlabel(r"$\log(M_*)$")
ax2b.set_xlim(11.0, 13.0)
ax2b.figure.canvas.draw()
# ax1b.set_title("z = 2")
ax1b.set_title("log(T/K) = 6")

lgd2 = legend.legend(ax1b)
lgd2.addLine((modellgd[0], 'black', lstyles[0], 1))
lgd2.addLine((modellgd[1], 'black', lstyles[1], 1))
lgd2.addLine((modellgd[2], 'black', lstyles[2], 1))
lgd2.loc = "upper right"
lgd2.fontsize = 12
lgd2.draw()

fig.subplots_adjust(left=0.12, top=0.8, right=0.95, wspace=0.05)

plt.savefig("/home/shuiyao_umass_edu/figures/tmp.pdf")

plt.show()

# plt.plot(msubs[::3], mstars[::3], "b.", alpha=0.4)
# plt.plot(x, y, "r-")


print ("DONE")
