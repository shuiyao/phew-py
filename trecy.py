import ioformat
import matplotlib.pyplot as plt
from scipy import sqrt, array, log10, linspace, histogram, isnan
import cosmology
import config_mpl
from matplotlib import gridspec
from pylab import setp
from matplotlib.mlab import griddata
import bin1d

DO_TMAX = False

flist = [\
         ["/scratch/shuiyao/scidata/newwind/p50n288fiducial/windsinfo.z2",\
          "/scratch/shuiyao/scidata/newwind/l25n144-phewoff/recyinfo.z2",\
          "/scratch/shuiyao/scidata/newwind/l25n144-phew/recyinfo.z2"],\
         ["/scratch/shuiyao/scidata/newwind/p50n288fiducial/windsinfo.z1",\
          "/scratch/shuiyao/scidata/newwind/l25n144-phewoff/recyinfo.z1",\
          "/scratch/shuiyao/scidata/newwind/l25n144-phew/recyinfo.z1"]\
        ]
clrs = ["magenta", "darkred", "green"]
lgds = ["Gadget3-L50", "Gadget3-L25", "GIZMO-L25"]
THUBBLE = log10(13.7e9)

def calc_trec(ai, af, log=True):
    for i in range(len(af)):
        if(af[i] > 1.0): af[i] = 1.0
    tf = cosmology.tcosmic(array(af))
    ti = cosmology.tcosmic(array(ai))
    trec = tf - ti
    if(log == True): trec = log10(trec)
    for i in range(len(af)):
        if(af[i] == 1.0): trec[i] = 20.0
    return trec

zmin = [2.0, 1.0]
MMIN, MMAX = 10.5, 13.0
NMBINS = 15
fig = plt.figure(1, figsize=(8,7))
gs = gridspec.GridSpec(3,2,height_ratios=[4,1,1])
patterns = []
for zi in [0, 1]:
    ax1 = plt.subplot(gs[zi])
    ax2 = plt.subplot(gs[zi+2])
    ax3 = plt.subplot(gs[zi+4])    
    for fi, fwind in enumerate(flist[zi]):
        if(DO_TMAX):
            fwindtmax = fwind[:-3]+"_tmax."+fwind[-2:]
            tmaxs = ioformat.rcol(fwindtmax, [0])
        x, y, xhot, xhotsub, yhot, yhotsub = [], [], [], [], [], []
        # ai, af, msi(FoF), msf(FoF), mstar(SKID), mvir(SO), flag(sopar)
        ai, af, mstar, mvir, flag = ioformat.rcol(fwind, [1, 2, 5, 6, 8], [4], linestart=1)
        zs = 1./array(ai) - 1.
        trec = calc_trec(ai, af, True)
        tlim = calc_trec([1./(zmin[zi]+1.)], [0.9999], True)
        for i in range(len(trec)):
            if(mvir[i] > MMIN and flag[i] == 1):
                # WARNING: In some earlier version, there's a bug that for ungrouped particles, mvir = -inf but mstar = mstar[-1]. So make sure use mvir[i] > MMIN is fine, but using mstar[i] would not be fine.
                x.append(mvir[i])
                y.append(trec[i])
                if(DO_TMAX):
                    if(tmaxs[i] > 5.5): # Winds that became hot
                        xhot.append(mvir[i])
                        yhot.append(trec[i])
                        if(trec[i] < tlim):
                            xhotsub.append(mvir[i])
                            yhotsub.append(trec[i])                        
        print "Number of winds: ", len(x)
        # x, y = bin1d.subsample(x, y, nonzero=True)
        xsub, ysub = [], []
        for i in range(len(x)):
            if(y[i] < tlim): # Stuff that recycles
                xsub.append(x[i])
                ysub.append(y[i])
        s = bin1d.bin1d(xsub, ysub, nbins=NMBINS, logbins=False, bounds=0.68)
        xline, yline = s.cen, s.median
        uline, lline = s.ubound, s.lbound
        idxs = (isnan(yline) == 0)
        xline, yline = list(array(xline)[idxs]), list(array(yline)[idxs])
        uline, lline = list(array(uline)[idxs]), list(array(lline)[idxs])
        ax1.plot(xline, yline, linestyle="-", color=clrs[fi])
        p, = ax1.fill(xline+xline[::-1], uline+lline[::-1], color=clrs[fi], alpha=0.4)
        if(DO_TMAX):
            s = bin1d.bin1d(xhotsub, yhotsub, nbins=NMBINS, logbins=False, bounds=0.68)
            xline, yline = s.cen, s.median
            uline, lline = s.ubound, s.lbound
            ax1.plot(xline, yline, linestyle=":", color=clrs[fi])        
        if(zi == 0):
            patterns.append(p)
        # ax1.plot(x[::10], y[::10], ".", color=clrs[fi], markersize=4, alpha=0.2)
        # -------- HISTOGRAM --------
        hist, edges = histogram(x, bins=linspace(MMIN, MMAX, NMBINS))
        hist_sub, edges = histogram(xsub, bins=linspace(MMIN, MMAX, NMBINS))
        if(DO_TMAX):
            hist_hot, edges = histogram(xhot, bins=linspace(MMIN, MMAX, NMBINS))
            hist_shot, edges = histogram(xhotsub, bins=linspace(MMIN, MMAX, NMBINS))
        # for j in range(len(hist))[1:]: hist[j] = float(hist[j] + hist[j-1])
        # hist = array(hist)/float(hist[-1])
        xmid, frac = [], []
        for i in range(len(hist)):
            if(hist[i] > 0):
                frac.append(float(hist_sub[i])/float(hist[i]))
                xmid.append(0.5*(edges[i] + edges[i+1]))
        ax2.plot(xmid, frac, color=clrs[fi], linestyle="-")

        if(DO_TMAX):
            fhot, fhotfrac = [], []
            for i in range(len(hist)):
                if(hist[i] > 0):
                    fhot.append(float(hist_hot[i])/float(hist[i]))
                    fhotfrac.append(float(hist_shot[i])/float(hist_hot[i]))
            ax2.plot(xmid, fhotfrac, color=clrs[fi], linestyle=":")        
            ax3.plot(xmid, fhot, color=clrs[fi], linestyle="-")        
        
    # ax = plt.gca()
    ax1.set_ylim(8.0, THUBBLE)
    # ax1.axhline(tlim[0], xmin=MMIN, xmax=MMAX, linestyle=":", color="black")
    ax1.plot([MMIN, MMAX], [tlim[0], tlim[0]], linestyle=":", color="black")
    ax1.text(12.4, tlim[0]+0.02, "z = %d" % (zmin[zi]), fontsize=12)
    ax2.set_ylim(0.0, 1.0)
    ax3.set_ylim(0.0, 1.0)        
    # ax1.set_ylim(8.5, tlim[0])
    # ax2.set_ylim(7.5, 10.25)    
    ax1.set_xlim(MMIN, MMAX)
    ax2.set_xlim(MMIN, MMAX)
    ax3.set_xlim(MMIN, MMAX)            
    # ax.set_xlabel("log(trec) [Gyr]")
    # ax.set_ylabel("Cumulative Fraction")
    ax3.set_xlabel(r"$\log(M_{vir}/M_\odot)$")
    ax1.yaxis.set_ticks([8.5, 9.0, 9.5, 10.0])    
    ax2.yaxis.set_ticks([0., 0.25, 0.5, 0.75, 1.0])
    ax3.yaxis.set_ticks([0., 0.25, 0.5, 0.75])
    ax3.xaxis.set_ticks([10.5, 11.0, 11.5, 12.0, 12.5])
    if(zi == 0):
        ax1.set_ylabel(r"$\log(t_{rec}) [yr]$")
        ax2.set_ylabel(r"$f_{rec}$")
        ax3.set_ylabel(r"$f_{hot}$")        
        ax1.legend(patterns, lgds, loc=3, fontsize=12)
    else:
        setp(ax1.get_yticklabels(),visible=False)
        setp(ax2.get_yticklabels(),visible=False)
        setp(ax3.get_yticklabels(),visible=False)                
    setp(ax1.get_xticklabels(),visible=False)
    setp(ax2.get_xticklabels(),visible=False)    
fig.subplots_adjust(wspace=0, hspace=0, left=0.15, right=0.9)
#plt.savefig("figures/trec.pdf")
plt.show()
