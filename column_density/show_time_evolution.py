import ioformat
from mymod import *
import matplotlib.pyplot as plt
from pltastro import frame, draw
from numpy import histogram, cumsum
import config_mpl

# x is perpendicular to the wind

fbase = "/scratch/shuiyao/Jneil/column_density/"
modelref = "T0.3_v1700_chi300_cond"
#modelref = "T0.3_v1700_chi300"
# modelref = "T0.3_v1000_chi300"

# model = "HC_v1000_chi300_cond"
# modelstr = "x300v1000"
# model = "T0.3_v1000_chi300"
# modelstr = "x300v1000kh"
# model = "T1_v1700_chi1000_cond"
# modelstr = "x1000v1700"
# model = "LowCond_v1700_chi300_cond"
# modelstr = "x300v1700c02"
model = "T0.3_v1700_chi300"
modelstr = "x300v1700kh"
mode = "frac" # frac or logN
fnamex = fbase + model + "_x.csv"
fnamey = fbase + model + ".csv"
# fname = fbase + "LowCond_v1700_chi300_cond_y.csv"
NPIX = 640000
# each pix: 2 pc x 2 pc
# Initial cloud size pi * 100 pc x 100 pc ~ 7854 pix
NPIXCLOUD = 7854

# rho = 1./2.89 * 1.0
rho = 2.0
mc = 1.0
Ncorr = log10(rho ** (2./3.) * mc ** (1./3.)) # Rho ** 2./3.
ncorr = rho ** (2./3.) * mc ** (-2./3.) # Rho ** 2./3.

clrs_t = ["magenta", "red", "orange", "yellow"]

def convert_list(lst):
    lst = array(lst)
    lst = lst[lst > 0]
    lst = log10(lst)
    return lst

def rcumsum(hist): # gadget3io/analyze/scratch.py
    if(mode == "frac"):
        cumfrac = cumsum(hist[::-1])[::-1]
        cumfrac = cumfrac / (float)(NPIXCLOUD)
    if(mode == "logN"):
        cumfrac = log10(hist)
    return cumfrac

# 1   2   3  4   5      6    7    8     9    10
# OVI CIV NV CII NeVIII CIII MgII SiIII SiIV HI

frm = frame.multi(4,2)
pars = frm.params
pars.figsize = (9, 9)
pars.left = 0.2
pars.top = 0.92
pars.bottom = 0.2
panels = frm.panels
panels.set_xlabels(r"$\log(N_{ion})$")
if(mode == "logN"):
    panels.set_ylims(0.0, 5.0)
    panels.set_yticks([0.0, 2.0, 4.0])
    panels.set_ylabels(r"$\log(n)$")
if(mode == "frac"):
    panels.set_ylims(0.0, 4.0)
    # panels.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])
    panels.set_yticks([0.0, 1.0, 2.0, 3.0, 4.0])    
    panels.set_ylabels(r"CS")
panels.set_xlims(10.0, 20.0)
panels.set_xticks([10.0, 12.0, 14.0, 16.0, 18.0])
fig, axs = draw(frm)

ionnames = ["XX","OVI","CIV","NV","CII","NeVIII","CIII","MgII","SiIII","SiIV","HI"]
ions = [1, 2, 7, 10]
nbins = linspace(10., 20., 101)
# Output Format: model.ionname
# Ncen t90x t75x t50x t25x t90y t75y t50y t25y

for axi, ion in enumerate(ions):
    cols = [ion, 10+ion, 20+ion, 30+ion]
    ionname = ionnames[ion]
    foutname = "histograms/"+model+"."+ionname
    foutnameref = "histograms/"+modelref+"."+ionname
    if(os.path.exists(foutname)):
        tab = genfromtxt(foutname, names=True)
        tabref = genfromtxt(foutnameref, names=True)        
        panel = 2 * axi
        axs[panel].plot(tab['Ncen'], rcumsum(tab['t90x']), "-", color=clrs_t[0])
        axs[panel].plot(tab['Ncen'], rcumsum(tab['t75x']), "-", color=clrs_t[1])
        axs[panel].plot(tab['Ncen'], rcumsum(tab['t50x']), "-", color=clrs_t[2])
        axs[panel].plot(tab['Ncen'], rcumsum(tab['t25x']), "-", color=clrs_t[3])
        axs[panel].plot(tabref['Ncen']+Ncorr, rcumsum(tabref['t90x'])/ncorr, "--", color=clrs_t[0])
        axs[panel].plot(tabref['Ncen']+Ncorr, rcumsum(tabref['t75x'])/ncorr, "--", color=clrs_t[1])
        axs[panel].plot(tabref['Ncen']+Ncorr, rcumsum(tabref['t50x'])/ncorr, "--", color=clrs_t[2])
        axs[panel].plot(tabref['Ncen']+Ncorr, rcumsum(tabref['t25x'])/ncorr, "--", color=clrs_t[3])        
        axs[panel].text(0.4, 0.85, ionname, fontsize=12, transform=axs[panel].transAxes)
        panel = 2 * axi + 1
        # axs[panel].plot(tab['Ncen']+Ncorr, rcumsum(tab['t50x']), "--", color=clrs_t[2])        
        axs[panel].plot(tab['Ncen'], rcumsum(tab['t90y']), "-", color=clrs_t[0])
        axs[panel].plot(tab['Ncen'], rcumsum(tab['t75y']), "-", color=clrs_t[1])
        axs[panel].plot(tab['Ncen'], rcumsum(tab['t50y']), "-", color=clrs_t[2])
        axs[panel].plot(tab['Ncen'], rcumsum(tab['t25y']), "-", color=clrs_t[3])
        axs[panel].plot(tabref['Ncen']+Ncorr, rcumsum(tabref['t90y'])/ncorr, "--", color=clrs_t[0])
        axs[panel].plot(tabref['Ncen']+Ncorr, rcumsum(tabref['t75y'])/ncorr, "--", color=clrs_t[1])
        axs[panel].plot(tabref['Ncen']+Ncorr, rcumsum(tabref['t50y'])/ncorr, "--", color=clrs_t[2])
        axs[panel].plot(tabref['Ncen']+Ncorr, rcumsum(tabref['t25y'])/ncorr, "--", color=clrs_t[3])        
        axs[panel].text(0.4, 0.85, ionname, fontsize=12, transform=axs[panel].transAxes)
    else:
        print "Creating New Files."
        fout = open(foutname, "w")
        fout.write("#Ncen t90x t75x t50x t25x t90y t75y t50y t25y\n")
        ioncds = ioformat.rcol(fnamex, cols, separator=",", linestart=1)
        hx = []
        for ioni in range(len(ioncds)): # t90 t75 t50 t25
            ioncds[ioni] = convert_list(ioncds[ioni])
            hist, edges = histogram(ioncds[ioni], bins=nbins)
            hx.append(hist)
        ioncds = ioformat.rcol(fnamey, cols, separator=",", linestart=1)
        hy = []
        for ioni in range(len(ioncds)): # t90 t75 t50 t25
            ioncds[ioni] = convert_list(ioncds[ioni])
            hist, edges = histogram(ioncds[ioni], bins=nbins)
            hy.append(hist)
        cen = 0.5 * (edges[1:] + edges[:-1])
        for i in range(len(cen)):
            fout.write("%5.2f %6d %6d %6d %6d %6d %6d %6d %6d\n" % \
                       (cen[i], hx[0][i], hx[1][i], hx[2][i], hx[3][i], \
                        hy[0][i], hy[1][i], hy[2][i], hy[3][i]))
        fout.close()

axs[0].set_title(modelstr+", x")
axs[1].set_title(modelstr+", y")
plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
plt.show()

print "done"
