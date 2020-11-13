import matplotlib.pyplot as plt
import ioformat
from scipy import logspace, array, linspace, sqrt, histogram, log10, histogram2d, ndimage, pi, median, meshgrid
import bin1d
from astroconst import pc, ac
from matplotlib import gridspec
from pylab import setp
from matplotlib.mlab import griddata
import matplotlib as mpl
import matplotlib.patches as mpatches
from matplotlib.colors import LogNorm
import config_mpl

FFORMAT = "NEW"
FBASE = "/scratch/shuiyao/scidata/newwind/"
simname = "l25n144-phewoff"
HALOMASS_CORRECTION = True
REDSHIFT = 2.0
if(REDSHIFT == 1.0):
    fwind = FBASE+simname+"/windsinfo.z1"
if(REDSHIFT == 2.0):
    fwind = FBASE+simname+"/windsinfo.z2"
XMIN, XMAX = 30., 600.
YMIN, YMAX = 50., 2000.
XBINS, YBINS = 20., 20.
CONTLEVELS = 5
SAVE_FIGURE = False
print "compiled."

def plotmedian(x, y, clr, opt=1, nbins=10, ax=[], xmin=30., xmax=150):
    xmid, ymid, ylist = [], [], []
    dx = (xmax - xmin) / nbins
    for i in range(nbins):
        xmid.append(xmin + (0.5+i)*dx)
        ylist.append([])
        ymid.append(0.)
    for i in range(len(x)):
        idx = int((x[i] - xmin) / dx)
        if(idx < 0):
            idx = 0
        if(idx > nbins-1):
            idx = nbins-1
        ylist[idx].append(y[i])
    for i in range(nbins):
        if((len(ylist[i])==0 or max(ylist[i])<0.1) and i>0):
            ymid[i] = ymid[i-1]
            xmid[i] = xmid[i-1]
        else:
            ymid[i] = median(ylist[i])
    if(opt==1):
        ax.plot(xmid, ymid, ".-", color=clr)
    if(opt==2):
        plt.plot(xmid, ymid, ".-", color=clr)
    print xmid, ymid

def V25Vc(fname=fwind):
    fig = plt.figure(1, figsize=(8,8))
    gs = gridspec.GridSpec(2,1,height_ratios=[4,1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    if(FFORMAT == "OLD"):
        mvir, rvir, vinit, v25 = ioformat.rcol(fname, [1,2,3,4], linestart=1)
    else:
        mvir, rvir, vinit, v25, rreturn = ioformat.rcol(fname, [1,2,3,4,5], linestart=1)
    if(HALOMASS_CORRECTION == True):
        mvir = array(mvir) - log10(8.0)
        rvir = array(rvir) / 2.0
    vc = sqrt(pc.G * 10.**array(mvir) * ac.msolar / (array(rvir) * ac.kpc)) / 1.e5
    print "Mvir Range: ", min(mvir), max(mvir)
    print "Vc Range: ", min(vc), max(vc)
    # ax1.plot(vc[::100], v25[::100], "b.")
    # ax1.plot(vc2, v252, ".", color="teal")
    if(FFORMAT == "NEW"):
        # Rreturn < 0.25 Rvir!
        for i in range(len(v25)):
            if(0.0 < rreturn[i] < 0.25*rvir[i]):
                v25[i] = -1
    plotmedian(vc, v25, "blue",nbins=15, ax=ax1, xmin=XMIN, xmax=XMAX)
    plotmedian(vc, vinit, "red",nbins=15, ax=ax1, xmin=XMIN, xmax=XMAX)
    xbins = logspace(log10(XMIN), log10(XMAX), XBINS)
    ybins = logspace(log10(YMIN), log10(YMAX), YBINS)
    xgrid, ygrid = meshgrid(xbins, ybins)
    z, edx, edy = histogram2d(vc, v25, bins=[xbins,ybins])
    z = z + 0.1
    z = z.T
    zf = ndimage.gaussian_filter(z, sigma=0.3, order=0)
    # cont = ax1.contour(xbins[1:], ybins[1:], zf, 10, colors="black")
    # cont = ax1.contourf(xbins[1:], ybins[1:], zf, CONTLEVELS, cmap=plt.cm.Purples, norm=LogNorm(vmin=z.min(), vmax=z.max()))
    # cont = ax1.contour(xbins[1:], ybins[1:], zf, CONTLEVELS, cmap=plt.cm.Reds, norm=LogNorm(vmin=z.min(), vmax=z.max()))    
    # cont = ax1.contourf(xbins[1:], ybins[1:], zf, 6, cmap=plt.cm.Purples, vmin=z.min(), vmax=z.max())
    # ax1.plot(vc[::200], v25[::200], "k.", markersize=2)
    ax1.pcolor(xbins[1:], ybins[1:], zf, cmap=plt.cm.Purples, norm=LogNorm(vmin=z.min(), vmax=z.max()))
    xline = linspace(XMIN, XMAX, 100)
    y50line = 0.854 * xline ** 1.12
    y95line = 1.85 * xline ** 1.10
    ax1.plot(xline, y50line, "k-")
    ax1.plot(xline, y95line, "k--")
    ax2.set_xlabel(r"$V_c [km/s]$")
    ax1.set_ylabel(r"$V_{25} [km/s]$")
    ax1.set_xlim(XMIN, XMAX)
    ax1.set_ylim(YMIN, YMAX)
    # ax1.set_ylim(9.,12.)
    ax1.set_title(simname+", Z = "+str(REDSHIFT))
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.legend([r"$V_{25}$", r"$V_{init}$"], fontsize=16)
    fig.subplots_adjust(hspace=0)
    setp(ax1.get_xticklabels(),visible=False)
    # Count the fraction that makes NOT to R25
    FRAC_STEP = 25.
    frac = [0.] * int(XMAX/FRAC_STEP)
    Ncount = [0.] * int(XMAX/FRAC_STEP)
    x = linspace(10., XMAX, int(XMAX/FRAC_STEP)) - 5.
    for i in range(len(vc)):
        idx = int(vc[i]/FRAC_STEP)
        if(idx > int(XMAX/FRAC_STEP)-1):
            idx = int(XMAX/FRAC_STEP)-1
        Ncount[idx] += 1.
        if(v25[i] != -1):
            frac[idx] += 1.
        # if(FFORMAT == "NEW"):
        #     if(0. < rreturn[i] < 0.25 * rvir[i]):
        #         frac[idx] -= 1. # Rreturn < 0.25 Rvir!
    Ntot = sum(Ncount)
    for i in range(len(frac)):
        if(Ncount[i] > 0.):
            frac[i] = frac[i] / Ncount[i]
    ax2.plot(x, frac, "k.-")
    # ax2.plot(x, array(Ncount)/Ntot, "-", color="teal")
    Ncount_Norm = 2.*max(Ncount)/Ntot
    ax2.bar(x, (array(Ncount)/Ntot)/Ncount_Norm, align="center", width=0.8*(x[1]-x[0]), color="grey")
    ax2.set_xlim(XMIN,XMAX)
    ax2.set_ylim(0.,1.1)
    ax2.set_xscale("log")
    ax2.xaxis.set_ticks([XMIN, 50., 100., 200., XMAX])
    ax2.xaxis.set_ticklabels([str(XMIN), "50", "100", "200", str(XMAX)])
    ax2.yaxis.set_ticks([0., 0.2, 0.4, 0.6, 0.8])
    ax2.set_ylabel(r"$f(R_{25} < R_{return})$")
    # plt.savefig("/home/shuiyao/161123V25VC.pdf")
    plt.show()
