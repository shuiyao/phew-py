# Make 2D map of when/where stars were formed.

from cosmology import acosmic, tcosmic
import ioformat
from numpy import genfromtxt, linspace, array, sort, insert
from numpy import where, inf, log10, isinf
import config_mpl
import matplotlib.pyplot as plt
import matplotlib as mpl

# WARNING: Line 49 _old !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Note about Tmax. Some cases Tmax < 0 for winds

print "compiled."
MODELNAME = "p50n288fiducial"
TITLENAME = "p50n288-fi"
def find_fnames(modelname, mstr):
    snapstr = "108"
    # mstr = "mh11"
    fname = "/scratch/shuiyao/scidata/gadget3io/"+modelname+"/"+modelname+"_"+snapstr+".starinfo."+mstr
    return fname
FIGNAME = "sfevents_"+MODELNAME+".pdf"

def arraysum(arr): # Turn histogram to cumulative distribution
    for i in range(len(arr))[1:]:
        arr[i] = arr[i] + arr[i-1]
    return arr

def read_starinfo(fname):
    stars = genfromtxt(fname, dtype='f8,f8,f8,f8,f8,f8,i8,i8', names=True)
    # a_form a_acc a_last Mass Mstar(at SF) Tmax GID HID
    # if Mvir < 0, it's a satellite galaxy
    return stars

def load_central_stars(fname):
    stars = read_starinfo(fname)
    nstars = len(stars)
    stars = stars[stars['GID'] == stars['HID']] # Only central galaxies
    print "Stars from central galaxies: %d/%d (%4.1f%%)" % \
        (len(stars), nstars, (float)(len(stars)) / (float)(nstars) * 100.)
    return stars

galsm11 = ("m11", 11.0, 11.5)
galsm12 = ("m12", 11.85, 12.15)
galsm13 = ("m13", 12.85, 13.15)
galss = [galsm11, galsm12, galsm13]
redz = ioformat.rcol("redshifts.txt", [0])
def show_lineage_of_selected_galaxies(modelname, gal_selection, ax=None):
    name, mmin, mmax = gal_selection[0], gal_selection[1], gal_selection[2]
    # folder = "/scratch/shuiyao/scidata/mergertree/"+modelname+"_old/"
    folder = "/scratch/shuiyao/scidata/mergertree/"+modelname+"/"    
    fin = open(folder + "lineages_" + name + ".dat", "r")
    xs = log10(array(redz) + 1)
    xarr, yarr = [], []
    for line in fin:
        spt = line.split()[1:]
        mstars = array(spt, dtype='float')[::-1]
        for i in range(len(mstars)):
            xarr.append(xs[i])
            yarr.append(mstars[i])
        # if(ax == None):
        #     plt.plot(xs, mstars, "r-", alpha=0.2)
        # else:
        #     ax.plot(xs, mstars, "r-", alpha=0.2)
    fin.close()
    plot1d.running_medians(xarr, yarr, ax=ax, nbins=50, bounds=0.68, scatter="bounds", lstyle=("red", '-', 1))

MMIN, MMAX = 8.0, 12.5
from pltastro import *
from pltastro import plot1d, plot2d
import pltastro
frm = pltastro.frame.multi(3, 3)
frm.params.height_ratios = [3, 3, 1]
frm.panels.ylabels[0] = r"$\log(M_*/M_\odot)$"
frm.panels.ylabels[3] = r"$z_{w}$"
frm.panels.ylabels[6] = r"$f(M_*)}$"
frm.panels.set_xlabels(r"$z_{sf}$")
frm.panels.set_ylims(MMIN, MMAX)
for i in [0,1,2]: # z - z
    frm.panels.yticks[i] = [9, 10, 11, 12]
for i in [3,4,5]: # z - z
    frm.panels.ylims[i] = (log10(7.0), 0.0)
    frm.panels.yticks[i] = [log10(5.0), log10(3.0), log10(2.0), 0.0]
    frm.panels.yticklabels[i] = ["4.0", "2.0", "1.0", "0.0"]
for i in [6,7,8]:
    frm.panels.ylims[i] = (0.0, 1.1)
    frm.panels.yticks[i] = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
# frm.panels.set_xticks([0.0, 0.33, 0.5, 0.8])
frm.panels.set_xticks([log10(5.0), log10(3.0), log10(2.0), 0.0])    
frm.panels.set_xticklabels(["4.0", "2.0", "1.0", "0.0"])
# frm.panels.set_xlims(0., 1.)
frm.panels.set_xlims(log10(7.0), 0.0)
frm.panels.set_xtickformat('%3.1f')
frm.params.figsize = (9, 9)
frm.params.bottom = 0.10
frm.params.title = TITLENAME
fig, axs = draw(frm)

mstrs = ['mh11', 'mh12', 'mh13']
titlestr = [r'$11.0 < \log(\frac{M_{vir}}{M_\odot}) < 11.5$',\
            r'$11.85 < \log(\frac{M_{vir}}{M_\odot}) < 12.15$',\
            r'$12.85 < \log(\frac{M_{vir}}{M_\odot}) < 13.15$']
lpats = [] # Legend patterns
# for fi in range(3):
for fi in [1]:
    mstr = mstrs[fi]
    fname = find_fnames(MODELNAME, mstr)
    # .starinfo.mh??
    stars = load_central_stars(fname)
    axs[fi].set_title(titlestr[fi])
    xarr, yarr, wt = stars['a_form'], stars['Mstar'], stars['Mass']
    masstot = sum(wt) # Total mass of stars in all a bins
    xarr = log10(1./xarr)
    percentlvls = [0.0, 0.10, 0.25, 0.5, 0.75]
    zlim = plot2d.plot2d(xarr, yarr, [0.0, 1.0], [MMIN, MMAX], ax=axs[fi], nbins=(40, 30), zarr=wt, method="contourf", logscale=False, clevels=percentlvls, percentiles=True)
    show_lineage_of_selected_galaxies(modelname=MODELNAME, gal_selection=galss[fi], ax=axs[fi])
    # xaxis is log(1+z)
    
    plot1d.histogram(xarr, wt, ax=axs[fi+6], nbins=30, xlim=[0.0, 1.0], lstyle=("black","-",1), cumulative=True, normalization=1.0, reverse=True)
    stars = stars[stars['a_last'] < 0] # Winds
    hotw = stars[abs(stars['Tmax']) > 5.5]
    coldw = stars[abs(stars['Tmax']) < 5.5]
    # Hot Wind Accretion
    xarr, yarr, y2arr, wt = log10(1./hotw['a_form']), hotw['Mstar'], log10(-1./hotw['a_last']), hotw['Mass']
    plot2d.plot2d(xarr, yarr, [0.0, 1.0], [MMIN, MMAX], ax=axs[fi], nbins=(40, 30), zarr=wt, method="contour", logscale=False, clevels=percentlvls, percentiles=True, lstyle=("yellowgreen", "-", 1), verbose=False)
    plot2d.plot2d(xarr, y2arr, [0.0, 1.0], [0.0, 1.0], ax=axs[fi+3], nbins=(40, 30), zarr=wt, method="contour", logscale=False, clevels=percentlvls, percentiles=True, lstyle=("yellowgreen", "-", 1), verbose=False)
    plot1d.histogram(xarr, wt, ax=axs[fi+6], nbins=30, xlim=[0.0, 1.0], lstyle=("yellowgreen","-",1), cumulative=True, normalization=sum(wt)/masstot, reverse=True)
    # Cold Wind Accretion
    xarr, yarr, y2arr, wt = log10(1./coldw['a_form']), coldw['Mstar'], log10(-1./coldw['a_last']), coldw['Mass']
    plot2d.plot2d(xarr, yarr, [0.0, 1.0], [MMIN, MMAX], ax=axs[fi], nbins=(40, 30), zarr=wt, method="contour", logscale=False, clevels=percentlvls, percentiles=True, lstyle=("darkgreen", "-", 1), verbose=False)
    plot2d.plot2d(xarr, y2arr, [0.0, 1.0], [0.0, 1.0], ax=axs[fi+3], nbins=(40, 30), zarr=wt, method="contour", logscale=False, clevels=percentlvls, percentiles=True, lstyle=("darkgreen", "-", 1), verbose=False)
    plot1d.histogram(xarr, wt, ax=axs[fi+6], nbins=30, xlim=[0.0, 1.0], lstyle=("darkgreen","-",1), cumulative=True, normalization=sum(wt)/masstot, reverse=True)
    
    axs[fi].plot([log10(3.), log10(3.)], [MMIN, MMAX], ":", color="lightgrey")
    axs[fi].plot([log10(2.), log10(2.)], [MMIN, MMAX], ":", color="lightgrey")
    axs[fi+3].plot([log10(3.), log10(3.)], [log10(7.0), 0.0], ":", color="lightgrey")
    axs[fi+3].plot([log10(2.), log10(2.)], [log10(7.0), 0.0], ":", color="lightgrey")
    axs[fi+3].plot([log10(7.0), 0.0], [log10(3.), log10(3.)], ":", color="lightgrey")
    axs[fi+3].plot([log10(7.0), 0.0], [log10(2.), log10(2.)], ":", color="lightgrey")
    axs[fi+6].plot([log10(3.), log10(3.)], [0.0, 1.1], ":", color="lightgrey")
    axs[fi+6].plot([log10(2.), log10(2.)], [0.0, 1.1], ":", color="lightgrey")        

    # axcbar = fig.add_axes([0.10+0.27*fi,0.06,0.24,0.015])
    # norm1 = mpl.colors.Normalize(vmin=zlim[0], vmax=zlim[1])
    # cdcbar = mpl.colorbar.ColorbarBase(axcbar, cmap=plt.cm.Purples, norm=norm1, orientation="horizontal")
    # cdcbar.set_ticks(zlim)
    # cdcbar.set_ticklabels(["0","max"])

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
lgds = [Patch(facecolor="purple", edgecolor='black', label="All"),\
        Line2D([0], [0], color="yellowgreen", linestyle="-", label="HotW"),\
        Line2D([0], [0], color="darkgreen", linestyle="-", label="ColdW")]
axs[0].legend(handles=lgds, loc="upper left")
# axs[0].legend(lgds, loc="upper left")
# plt.savefig(FIGNAME)
plt.show()
