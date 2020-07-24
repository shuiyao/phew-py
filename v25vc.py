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
simname = "p50n288ezwc"
#simname = "p50n288zw"
#simname = "p50n288sw"
REDSHIFT = 1.0
# XMIN, XMAX = 30., 400.
# YMIN, YMAX = 50., 2000.
XMIN, XMAX = 30., 800.
YMIN, YMAX = 50., 2500.
# XMIN, XMAX = 30., 800. # Whole range
# YMIN, YMAX = 50., 3000.
XBINS, YBINS = 20., 20.
CONTLEVELS = 5
print "compiled."
SAVEFIG = False

def get_fname(simname, z=REDSHIFT):
    if(z == 1.0):
        fwind = FBASE+simname+"/windsinfo.z1"
    if(z == 2.0):
        fwind = FBASE+simname+"/windsinfo.z2"
    return fwind

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

def V25Vc(fname, ax1, ax2, titlename):
    import cosmology
    if(FFORMAT == "OLD"):
        mvir, rvir, vinit, v25 = ioformat.rcol(fname, [1,2,3,4], linestart=1)
    else:
        mvir, rvir, vinit, v25, rreturn = ioformat.rcol(fname, [1,2,3,4,5], linestart=1)    
    # vc = sqrt(pc.G * 10.**array(mvir) * ac.msolar / (array(rvir) * ac.kpc)) / 1.e5
    vc = []
    for i in range(len(mvir)):
        vc.append(cosmology.Vc(10.**mvir[i]*ac.msolar, REDSHIFT) / 1.e5)
    print "Mvir Range: ", min(mvir), max(mvir)
    print "Vc Range: ", min(vc), max(vc)
    # ax1.plot(vc[::100], v25[::100], "b.")
    # ax1.plot(vc2, v252, ".", color="teal")
    if(FFORMAT == "NEW"):
        # Rreturn < 0.25 Rvir!
        for i in range(len(v25)):
            if(0.0 < rreturn[i] < 0.25*rvir[i]):
                v25[i] = -1
    plotmedian(vc, v25, "cyan",nbins=15, ax=ax1, xmin=XMIN, xmax=XMAX)
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
    xsig16 = array([50., 85.])
    ysig16 = 0.7 * xsig16 ** 1.60
    ax1.plot(xline, y50line, "k-")
    ax1.plot(xline, y95line, "k--")
    ax1.plot(xsig16, ysig16, "g-")    
    ax2.set_xlabel(r"$V_c [km/s]$")
    ax1.set_xlim(XMIN, XMAX)
    ax1.set_ylim(YMIN, YMAX)
    # ax1.set_ylim(9.,12.)
    ax1.set_title(titlename+", z = "+str(REDSHIFT))
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
    print "Ncount = ", Ncount
    print "frac = ", frac
    # ax2.plot(x, array(Ncount)/Ntot, "-", color="teal")
    Ncount_Norm = 2.*max(Ncount)/Ntot
    ax2.bar(x, (array(Ncount)/Ntot)/Ncount_Norm, align="center", width=0.8*(x[1]-x[0]), color="grey")
    ax2.set_xlim(XMIN,XMAX)
    ax2.set_ylim(0.,1.1)
    ax2.set_xscale("log")
    ax2.yaxis.set_ticks([0., 0.2, 0.4, 0.6, 0.8, 1.0])
    ax2.axhline(1.0, XMIN, XMAX, color="grey", linestyle=":")

simnames = ["p50n288zw", "l50n288-phewoff", "l25n288-phewoff-fw"]
titles = ["Gadget3", "GIZMO-L50", "GIZMO-L25"]
# simnames = ["l50n288-phewoff", "l50n288-phew-m5", "l25n288-phew-m5"]
# titles = ["GIZMO", "PhEW-Mc5-L50", "PhEW-Mc5-L25"]
# simnames = ["l50n288-phewoff", "l50n288-phew-m4", "l25n288-phew-m4"]
# titles = ["GIZMO", "PhEW-Mc4-L50", "PhEW-Mc4-L25"]
fig = plt.figure(1, figsize=(10,8))
gs = gridspec.GridSpec(2,3,height_ratios=[4,1])
for si, simname in enumerate(simnames):
    ax1 = plt.subplot(gs[si])
    ax2 = plt.subplot(gs[si+3])
    V25Vc(fname=get_fname(simname), ax1=ax1, ax2=ax2, titlename=titles[si])
    if(si == 0): # Leftmost panel
        ax1.set_ylabel(r"$V_{25} [km/s]$")
        ax2.set_ylabel(r"$f(R_{25} < R_{return})$")
        ax2.xaxis.set_ticks([XMIN, 50., 100., 200., XMAX])        
        ax2.xaxis.set_ticklabels([str(XMIN), "50", "100", "200", str(XMAX)])
    else:
        setp(ax1.get_yticklabels(),visible=False)
        setp(ax2.get_yticklabels(),visible=False)
        ax2.xaxis.set_ticks([50., 100., 200., XMAX])
        ax2.xaxis.set_ticklabels(["50", "100", "200", str(XMAX)])    

fig.subplots_adjust(wspace=0, left=0.1, right=0.9)
if(SAVEFIG == True):
    if(REDSHIFT == 1.0):
        plt.savefig("./figures_final/v25vc_z1.pdf")
    if(REDSHIFT == 2.0):
        plt.savefig("./figures_final/v25vc_z2.pdf")
plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
plt.show()
