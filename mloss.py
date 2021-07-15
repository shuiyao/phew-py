import matplotlib.pyplot as plt
import ioformat
from scipy import logspace, array, linspace, sqrt, histogram, log10, histogram2d, ndimage, pi, median, meshgrid
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
from myinit import *

# USEFUL FUNCTIONS:
# plot_mloss_radius()

# PARENT: /scratch/shuiyao/sci/newwind/mloss/mloss.py

modelname = "l50n576-phew-m5"
MC_INIT = 2.0e38
MMIN, MMAX = 11.0, 13.5
#MMIN, MMAX = 11.0, 12.5
#unit_m = 433697.735404
unit_m = 433697.735404 * 8.0

REDSHIFT = 1.0
if(REDSHIFT == 2.0):
    subwfolder = "z2"
if(REDSHIFT == 1.0):
    subwfolder = "z1"

ascale = 1./(REDSHIFT+1.)
# Mgasp = 82656250.0 # P6N36
Mgasp = 9.3e7
outputbase = DIRS['SCIDATA']+"phew/"
fwind = outputbase+modelname+"/"+"phewsinfo."+subwfolder
fmloss = outputbase+modelname+"/"+"mlossinfo."+subwfolder
alphavalue = 0.4
boundsvalue = 0.6

mclast_cut = 0.1
print("compiled")

# phewp = genfromtxt(fmloss, names=True)
# modelname = "l50n288-phew-m5"
# fmloss2 = outputbase+modelname+"/"+"mlossinfo."+subwfolder
# phewp2 = genfromtxt(fmloss2, names=True)
# plt.plot(phewp['Mvir'], phewp['Rlast'], "b.", alpha=0.2)
# plt.plot(phewp2['Mvir'], phewp2['Rlast'], "r.", alpha=0.2)
# plt.show()

def Mclast():
    fig = plt.figure(1, figsize=(6,9))
    gs = gridspec.GridSpec(2,1,height_ratios=[4,1])
    ax = plt.subplot(gs[0])
    mvir, mclast, mratio = ioformat.rcol(fwind, [1,5,14], linestart=1)
    print("Mvir Range: ", min(mvir), max(mvir))
    x, y = [], []
    for i in range(len(mclast)):
        if(mclast[i] > mclast_cut and mratio[i] < 20.0):
            x.append(mvir[i])
            y.append(mclast[i])
    xbins = linspace(11.0, 13.5, 30)
    ybins = linspace(0.05, 1.0, 40)
    xgrid, ygrid = meshgrid(xbins, ybins)
    z, edx, edy = histogram2d(x, y, bins=[xbins,ybins])
    z = z.T + 0.01
    zf = ndimage.gaussian_filter(z, sigma=1.0, order=0)
    cont = ax.contour(xbins[1:], ybins[1:], zf, colors="red")
    #plt.pcolor(xgrid, ygrid, z, cmap="Purples", norm=LogNorm(vmin=z.min(), vmax=z.max()))
    ax.pcolor(xgrid, ygrid, z, cmap="Purples")
    setp(ax.get_xticklabels(), visible=False)
    ax.set_ylabel("Mc (Rejoin)")
    plt.title(modelname+", Z~1.0")
    ax = plt.subplot(gs[1])
    plt.subplots_adjust(hspace=0.0, top=0.9, bottom=0.15)
    hist1, bins = histogram(mvir, bins=linspace(11.0,13.5,30))
    hist2, bins = histogram(x, bins=linspace(11.0,13.5,30))
    hist = []
    for i in range(len(hist1)):
        if(hist1[i] > 0): hist.append(float(hist2[i])/float(hist1[i]))
        else: hist.append(0.0)
    #width=0.8*(bins[1]-bins[0])
    center = (bins[:-1] + bins[1:]) / 2.0
    ax.plot(center, hist, "b.-")
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("Mvir")
    ax.set_ylabel("f_rej")
    plt.show()

def plot_last_signal():
    mvir, rvir, rlast, mclast = ioformat.rcol(fwind, [8, 10, 6, 5], linestart=1)
    mbin_min = [11.25, 11.75, 12.25, 12.75]
    mbin_max = [11.75, 12.25, 12.75, 13.75]
    clrs = ["magenta", "red", "orange", "yellow"]
    wpmbin = []
    class mvir_bin():
        def __init__(self, mmin, mmax):
            self.mmin = mmin
            self.mmax = mmax
            self.rlast = []
            self.rratio = []
            self.mclast = []

    ax = plt.figure(1, figsize=(6,6)).add_subplot(111)        
    for i in range(len(mbin_min)):
        wpmbin.append(mvir_bin(mbin_min[i], mbin_max[i]))
        for j in range(len(mvir)):
            if (mbin_min[i] < mvir[j] < mbin_max[i]):
                wpmbin[-1].rlast.append(rlast[j])
                wpmbin[-1].mclast.append(mclast[j])
                wpmbin[-1].rratio.append(rlast[j]/rvir[j])
        # x, y = bin1d.subsample(wpmbin[-1].mclast, wpmbin[-1].rlast, nonzero=True)
        # print len(wpmbin[-1].mclast), len(wpmbin[-1].rratio)
        if(len(wpmbin[-1].mclast) > 0):
            x, y = bin1d.subsample(wpmbin[-1].mclast, wpmbin[-1].rratio, nonzero=True)
            s = bin1d.bin1d(x, y, nbins=10, bounds=boundsvalue)
            xline, yline = s.value, s.median
            uline, lline = s.ubound, s.lbound
            p, = ax.plot(xline, yline, "-", color=clrs[i]) # Draw line for all models    
            pf, = plt.fill(xline+xline[::-1], uline+lline[::-1], color=clrs[i], alpha=alphavalue)
    plt.ylabel(r'$r_{last} / R_{vir}$')
    plt.xlabel(r'$Mc_{last}$')
    legends = []
    legends.append(mpatches.Patch(color="yellow"))
    legends.append(mpatches.Patch(color="orange"))
    legends.append(mpatches.Patch(color="red"))
    legends.append(mpatches.Patch(color="magenta"))
    plt.legend(legends, [r'$12.75<M_{vir}<13.75$',r'$12.25<M_{vir}<12.75$',r'$11.75<M_{vir}<12.25$',r'$11.25<M_{vir}<11.75$'], loc=1, fontsize=12)
    plt.show()

def plot_mloss_radius():
    fig = plt.figure(1, figsize=(7,9))
    gs = gridspec.GridSpec(3,1,height_ratios=[1,1,1])
    ax0 = plt.subplot(gs[0])
    Massloss_time(ax=ax0)
    plt.title(modelname+", z~"+str(REDSHIFT)[:3])
    ax1 = plt.subplot(gs[1])
    Massloss_radius(ax=ax1, rnorm=False)
    ax2 = plt.subplot(gs[2])
    Massloss_radius(ax=ax2, rnorm=True)
    setp(ax0.get_xticklabels(), visible=False)
    setp(ax1.get_xticklabels(), visible=False)
    fig.subplots_adjust(hspace=0, left=0.15)
    plt.savefig(DIRS['FIGURE']+'tmp.pdf')
    plt.show()

def plot_mloss_velocity():
    fig = plt.figure(1, figsize=(6,9))
    gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
    ax1 = plt.subplot(gs[0])
    Massloss_velocity(ax=ax1, rnorm=False)
    plt.title(modelname+", z~1.0")
    ax2 = plt.subplot(gs[1])
    Massloss_velocity(ax=ax2, rnorm=True)
    setp(ax0.get_xticklabels(), visible=False)
    setp(ax1.get_xticklabels(), visible=False)    
    fig.subplots_adjust(hspace=0)
    plt.show()

def Massloss_radius(fin=fmloss, rnorm=False, ax=[]):
    mclast_cut = 0.1
    phewp = genfromtxt(fin, names=True)
    # SELECT ONLY THOSE ANNIHILATED
    phewp = phewp[phewp['Rlast'] > 0]
    r75, r50, r25, rlast = phewp['R75'], phewp['R50'], phewp['R25'], phewp['Rlast']
    mvir = log10(phewp['Mvir'])
    if(rnorm==True):
        r25 /= phewp['Rvir']
        r50 /= phewp['Rvir']
        r75 /= phewp['Rvir']
        rlast /= phewp['Rvir']
    print(max(mvir), min(mvir))

    NBIN = 10
    # ax.plot(mvir[::10], rlast[::10], "k.", markersize=6)
    plotmedian(mvir, rlast, ax, nbins=NBIN, clr="yellow", alphavalue=alphavalue)    
    plotmedian(mvir, r25, ax, nbins=NBIN, clr="orange", alphavalue=alphavalue, verbose=True)
    plotmedian(mvir, r50, ax, nbins=NBIN, clr="red", alphavalue=alphavalue)
    plotmedian(mvir, r75, ax, nbins=NBIN, clr="magenta", alphavalue=alphavalue)

    ax.set_xlim(MMIN, MMAX)
    ax.set_xlabel(r'$Log(M_{vir}/M_\odot)$')
    if(rnorm == True):
        ax.set_ylabel(r'$r/R_{vir}$')
        ax.set_ylim(0.0,1.5)
        ax.set_yticks([0., 0.5, 1.0])
        # ax.plot([MMIN, MMAX], [0.25, 0.25], "k:")
    else:
        ax.set_ylabel("r [kpc]")
        ax.set_ylim(0.,450.)
        ax.set_yticks([0., 100., 200., 300., 400.])        
    legends = []
    legends.append(mpatches.Patch(color="yellow"))
    legends.append(mpatches.Patch(color="orange"))
    legends.append(mpatches.Patch(color="red"))
    legends.append(mpatches.Patch(color="magenta"))
    plt.legend(legends, [r'$r_{last}$', r'$r_{25}$', r'$r_{50}$', r'$r_{75}$'], loc="upper right", fontsize=16)
    # ax.text(0.8, 0.11, "+ V10/Vinit", color="blue", transform=ax.transAxes)

def Massloss_time(fin=fmloss, ax=[]):
    mclast_cut = 0.1
    phewp = genfromtxt(fin, names=True)
    # SELECT ONLY THOSE ANNIHILATED
    phewp = phewp[phewp['Rlast'] > 0]
    t75, t50, t25, tlast = phewp['t75'], phewp['t50'], phewp['t25'], phewp['tlast']
    mvir = log10(phewp['Mvir'])

    NBIN = 10
    plotmedian(mvir, tlast, ax, nbins=NBIN, clr="yellow", alphavalue=alphavalue)    
    plotmedian(mvir, t25, ax, nbins=NBIN, clr="orange", alphavalue=alphavalue)
    plotmedian(mvir, t50, ax, nbins=NBIN, clr="red", alphavalue=alphavalue)
    plotmedian(mvir, t75, ax, nbins=NBIN, clr="magenta", alphavalue=alphavalue)

    ax.set_xlim(MMIN, MMAX)
    ax.set_xlabel(r'$Log(M_{vir}/M_\odot)$')
    ax.set_ylabel(r'$Time [Myr]$')
    ax.set_ylim(0.,1100)
    ax.set_yticks([0., 250., 500., 750., 1000.])
    legends = []
    legends.append(mpatches.Patch(color="yellow"))
    legends.append(mpatches.Patch(color="orange"))
    legends.append(mpatches.Patch(color="red"))
    legends.append(mpatches.Patch(color="magenta"))
    plt.legend(legends, [r'$t_{last}$', r'$t_{25}$', r'$t_{50}$', r'$t_{75}$'], loc="upper right", fontsize=16)
    # ax.text(0.8, 0.11, "+ V10/Vinit", color="blue", transform=ax.transAxes)

def Massloss_velocity(fin=fwind, rnorm=False, ax=[]):
    mclast_cut = 0.1
    # ax = plt.figure(1, figsize=(8,8)).add_subplot(111)
    # mvir, vinit, v25, v50, v75, rlast, mclast = ioformat.rcol(fin, [1,2,7,8,9,6,5], linestart=1)
    mvir, rvir, vinit = [], [], []
    vc, v25, v50, v75, vlast = [], [], [], [], []
    f = open(fin, "r")
    f.readline()
    for line in f:
        spt = line.split()
        if(float(spt[5]) <= mclast_cut and float(spt[14]) < 15.0):
            mvir.append(float(spt[1]))
            rvir.append(float(spt[2]))            
            vinit.append(float(spt[3]))
            v25.append(float(spt[10]))
            v50.append(float(spt[11]))
            v75.append(float(spt[12]))
            vlast.append(float(spt[13]))
            vc.append(sqrt(pc.G*10.**mvir[-1]*ac.msolar/(rvir[-1]*ac.kpc)) / 1.e5)
    if(rnorm==True):
        vc = array(vc)
        v25 /= vc
        v50 /= vc
        v75 /= vc
        vlast /= vc
        print(max(mvir), min(mvir))

    plotmedian(mvir, vlast, ax, nbins=20, clr="yellow", alphavalue=alphavalue)    
    plotmedian(mvir, v75, ax, nbins=20, clr="orange", alphavalue=alphavalue)
    plotmedian(mvir, v50, ax, nbins=20, clr="red", alphavalue=alphavalue)
    plotmedian(mvir, v25, ax, nbins=20, clr="magenta", alphavalue=alphavalue)

    ax.set_xlim(11.5,13.5)
    ax.set_xlabel(r'$Log(M_{vir}/M_\odot)$')
    if(rnorm == True):
#        ax.set_ylabel(r'$v/V_{init}$')
        ax.set_ylabel(r'$v/V_{c}$')        
        ax.set_ylim(0.,3.5)        
    else:
        ax.set_ylabel("v [km/s]")
        ax.set_ylim(80.,1500.)
        ax.set_yscale("log")        
    legends = []
    legends.append(mpatches.Patch(color="yellow"))    
    legends.append(mpatches.Patch(color="orange"))
    legends.append(mpatches.Patch(color="red"))
    legends.append(mpatches.Patch(color="magenta"))
    plt.legend(legends, [r'$v_{last}$', r'$v_{75}$', r'$v_{50}$', r'$v_{25}$'], loc=4)

def plotmedian(x0, y0, ax, nbins=20, clr="blue", alphavalue=0.4, verbose=False):
    x, y = bin1d.subsample(x0, y0, nonzero=True)
    s = bin1d.bin1d(x, y, nbins=nbins, bounds=boundsvalue)
    xline, yline = s.cen, s.median
    uline, lline = s.ubound, s.lbound
    if(verbose == True):
        print(xline, yline)
    for i in range(len(yline))[1:-1]:
        if(yline[i] == 0):
            yline[i] = 0.5 * (yline[i-1] + yline[i+1])
            uline[i] = 0.5 * (uline[i-1] + uline[i+1])
            lline[i] = 0.5 * (lline[i-1] + lline[i+1])            
    p, = ax.plot(xline, yline, "-", color=clr) # Draw line for all models    
    pf, = plt.fill(xline+xline[::-1], uline+lline[::-1], color=clr, alpha=alphavalue)
    
