# A specific plot for Neal.
# The r/rvir of PhEWs at 75% mass loss

import matplotlib.pyplot as plt
import ioformat
from scipy import logspace, array, linspace, sqrt, histogram, log10, pi
import bin1d
from numpy import genfromtxt
from matplotlib.colors import LogNorm
from astroconst import pc, ac
from matplotlib import gridspec
from pylab import setp
from scipy import rand
import matplotlib as mpl
import matplotlib.patches as mpatches
import config_mpl

REDSHIFT = 1.0
if(REDSHIFT == 2.0):
    subwfolder = "z2"
if(REDSHIFT == 1.0):
    subwfolder = "z1"

ascale = 1./(REDSHIFT+1.)

models = ["l25n288-phew-m4","l25n288-phew-m5","l50n288-phew-m4","l50n288-phew-m5"]
MC_INIT = [2.0e37, 2.0e38, 2.0e37, 2.0e38]
unit_m = [433697.735404, 433697.735404, 433697.735404*8.0, 433697.735404*8.0]
clrs = ["forestgreen", "red", "lime", "orange"]
lgds = [r"$L25, M_c=10^4M_\odot$", r"$L25, M_c=10^5M_\odot$", r"$L50, M_c=10^4M_\odot$", r"$L50, M_c=10^5M_\odot$"]
alphavalue = 0.4
boundsvalue = 0.60

MMIN, MMAX = 11.0, 13.5
Mgasp = 9.3e7
outputbase = "/home/shuiyao_umass_edu/scidata/"

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

fig, ax = plt.subplots(1,1, figsize=(7,6))
legends = []

for mi, modelname in enumerate(models):
    fwind = outputbase+modelname+"/"+"phewsinfo."+subwfolder
    fmloss = outputbase+modelname+"/"+"mlossinfo."+subwfolder
    plt.title("z ~"+str(REDSHIFT)[:3])
    phewp = genfromtxt(fmloss, names=True)
    phewp = phewp[phewp['Rlast'] > 0]    
    r25 = phewp['R25'] / phewp['Rvir']
    mvir = phewp['Mvir']
    
    NBIN = 10
    plotmedian(mvir, r25, ax, nbins=NBIN, clr=clrs[mi], alphavalue=alphavalue, verbose=True)

    ax.set_xlim(MMIN, MMAX)
    ax.set_xlabel(r'$Log(M_{vir}/M_\odot)$')
    ax.set_ylabel(r'$r_{25}/R_{vir}$')
    ax.set_ylim(0.,1.5)
    ax.set_yticks([0., 0.5, 1.0, 1.5])

    legends.append(mpatches.Patch(color=clrs[mi]))

plt.legend(legends, lgds, loc="upper right")
fig.subplots_adjust(hspace=0, left=0.15)
plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
plt.show()

