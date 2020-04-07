import ioformat
import matplotlib.pyplot as plt
from scipy import sqrt, array, log10, linspace
from astropy.table import Table
from astropy.io import ascii
import cosmology
from matplotlib import gridspec
from pylab import setp

import config_mpl

bs16 = Table.read("models_bs16.dat", format="ascii", data_start=0)
bs16pred = Table.read("models_bs16_predictions.dat", format="ascii", data_start=0)

modelnames = ["models_191110", "models_191110_q95", "models_191110_q80"]
# modelnames = ["models_191110", "models_191110_nokhi", "models_191110_q80"]

FIGNAME = "./figures/mloss_qs.pdf"

def read_models(midx):
    ana = Table.read(modelnames[midx]+".dat", format="ascii", data_start=0)
    ana['t75']= ana['t75'] + 1
    ana['t50']= ana['t50'] + 1
    ana['t25']= ana['t25'] + 1
    return ana

symlst = ["*", "^", ".", ".", ".", "."]
# clrs = ["blue", "purple", "red", "green", "yellow", "purple", "brown"]
clrs = ["blue", "purple", "cyan"]
vmax = [150., 120., 250., 80., 210., 320.,]
dvlbl = [30., 20., 50., 15., 40., 60.,]

fig = plt.figure(1, figsize=(10,7))
# ax1 = fig.add_subplot(111)
# fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True)
# fig.subplots_adjust(hspace=0, wspace=0)
# setp(ax1.get_xticklabels(),visible=False)
# fig.suptitle("Title")
NPANELS = 8
gs = gridspec.GridSpec(NPANELS/2,2)
axs = []
for i in range(NPANELS):
    axs.append(plt.subplot(gs[i]))
    setp(axs[-1].get_xticklabels(),visible=False)
    setp(axs[-1].get_yticklabels(),visible=False)    
    axs[-1].set_xlim(0., 35.)
    axs[-1].set_ylim(0., 1.)    
    axs[-1].yaxis.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8])
    axs[-1].xaxis.set_ticks([0, 10, 20, 30])
    if((float)(i)/2.0 == (int)(i/2)):
        setp(axs[i].get_yticklabels(),visible=True)
        axs[i].set_ylabel(r"$\frac{M_c(t)}{M_c(t=0)}$")
setp(axs[-2].get_xticklabels(),visible=True)    
setp(axs[-1].get_xticklabels(),visible=True)
axs[-2].set_xlabel(r'$t/\tau_{cc}$')
axs[-1].set_xlabel(r'$t/\tau_{cc}$')
fig.subplots_adjust(hspace=0, wspace=0, top=0.9, bottom=0.1)

for midx in range(len(modelnames)):
    ana = read_models(midx)
    for i in range(len(ana)):
        if(i > 7): continue
        idx = ana['ID'][i]
        ms = 9
        axs[idx].plot(ana['t75'][i], 0.75, color=clrs[midx], marker=symlst[0], markersize=ms)
        axs[idx].plot(ana['t50'][i], 0.50, color=clrs[midx], marker=symlst[1], markersize=ms)
        axs[idx].plot(ana['t25'][i], 0.25, color=clrs[midx], marker=symlst[2], markersize=ms)
        x = [0.0, ana['t75'][i], ana['t50'][i], ana['t25'][i]]
        y = [1.0, 0.75, 0.50, 0.25]
        axs[idx].plot(x, y, "-", color=clrs[midx])
        # axs[idx].plot([0., 40.], [ana['vs'][i], ana['vs'][i]], "k--")

for i in range(NPANELS):
    ms = 9
    idx = bs16['ID'][i]
    if(idx > 7): continue
    axs[idx].plot(bs16['t75'][i], 0.75, color="black", marker=symlst[0], markersize=ms)
    axs[idx].plot(bs16['t50'][i], 0.50, color="black", marker=symlst[1], markersize=ms)
    axs[idx].plot(bs16['t25'][i], 0.25, color="black", marker=symlst[2], markersize=ms)
    x = [0.0, bs16['t75'][i], bs16['t50'][i], bs16['t25'][i]]
    y = [1.0, 0.75, 0.50, 0.25]
    axs[idx].plot(x, y, "-", color="black")

for i in range(NPANELS):
    idx = bs16pred['ID'][i]    
    x = [0.0, bs16pred['tevap'][i]]
    y = [1.0, 0.0]
    axs[idx].plot(x, y, "--", color="grey")

for i in range(NPANELS):
    axs[i].text(0.7, 0.1, bs16['Modelname'][i], color='black', fontsize=12, transform=axs[i].transAxes)

axs[0].plot(0.05, 0.40, color="black", marker=symlst[0], markersize=12, transform=axs[0].transAxes)
axs[0].plot(0.05, 0.25, color="black", marker=symlst[1], markersize=12, transform=axs[0].transAxes)
axs[0].plot(0.05, 0.10, color="black", marker=symlst[2], markersize=12, transform=axs[0].transAxes)
axs[0].text(0.08, 0.40, r'$t_{75}$', color='black', fontsize=12, transform=axs[0].transAxes, va='center')
axs[0].text(0.08, 0.25, r'$t_{50}$', color='black', fontsize=12, transform=axs[0].transAxes, va='center')
axs[0].text(0.08, 0.10, r'$t_{25}$', color='black', fontsize=12, transform=axs[0].transAxes, va='center')    

# xoff = 0.05
# for i in range(NPANELS):
#     axs[i].plot(0.80-xoff, 0.9, color="grey", marker=symlst[1], markersize=12, transform=axs[i].transAxes)
#     axs[i].text(0.83-xoff, 0.9, "BS16", color='grey', fontsize=12, transform=axs[i].transAxes, va='center')
#     axs[i].plot(0.80-xoff, 0.72, color=clrs[0], marker=symlst[1], markersize=12, transform=axs[i].transAxes)   
#     axs[i].text(0.83-xoff, 0.72, r'$\hat{q}=0.90$', color=clrs[0], fontsize=12, transform=axs[i].transAxes, va='center')
#     axs[i].plot(0.80-xoff, 0.54, color=clrs[1], marker=symlst[1], markersize=12, transform=axs[i].transAxes)   
#     axs[i].text(0.83-xoff, 0.54, r'$\hat{q}=0.95$', color=clrs[1], fontsize=12, transform=axs[i].transAxes, va='center')
#     axs[i].plot(0.80-xoff, 0.36, color=clrs[2], marker=symlst[1], markersize=12, transform=axs[i].transAxes)   
#     axs[i].text(0.83-xoff, 0.36, r'$\hat{q}=0.80$', color=clrs[2], fontsize=12, transform=axs[i].transAxes, va='center')

from pltastro import legend
lgd = legend.legend(axs[6])
lgd.addLine(("BS16 simulation", "black", "-", 1))
lgd.addLine(("BS16 prediction", "grey", "--", 1))
lgd.loc = "upper right"
lgd.fontsize = 11
lgd.draw()
lgd = legend.legend(axs[7])
lgd.addLine(("PhEW $\hat{q}=0.90$", clrs[0], "-", 1))
lgd.addLine(("PhEW $\hat{q}=0.95$", clrs[1], "-", 1))
lgd.addLine(("PhEW $\hat{q}=0.80$", clrs[2], "-", 1))
lgd.loc = "upper right"
lgd.fontsize = 11
lgd.draw()


# plt.savefig("./figures/mloss_nokhi.pdf")
plt.savefig(FIGNAME)    
plt.show()


