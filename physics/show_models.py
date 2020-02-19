import ioformat
import matplotlib.pyplot as plt
from scipy import sqrt, array, log10, linspace, pi
from astropy.table import Table
from astropy.io import ascii
import cosmology
from matplotlib import gridspec
import matplotlib.ticker as ticker
from pylab import setp

import config_mpl

SHOW_BS16PRED = False
bs16 = Table.read("models_bs16.dat", format="ascii", data_start=0)
bs16pred = Table.read("models_bs16_predictions.dat", format="ascii", data_start=0)
ana = Table.read("models_191110.dat", format="ascii", data_start=0)
ana_sph = Table.read("models_spherical_190912.dat", format="ascii", data_start=0)
ana['t75']= ana['t75'] + 1
ana['t50']= ana['t50'] + 1
ana['t25']= ana['t25'] + 1
ana_sph['t75']= ana_sph['t75'] + 1
ana_sph['t50']= ana_sph['t50'] + 1
ana_sph['t25']= ana_sph['t25'] + 1

symlst = ["*", "^", ".", ".", ".", "."]
clrs = ["blue", "purple", "red", "green", "yellow", "purple", "brown"]
vmax = [160., 120., 60., 60., 280., 280., 60., 60.]
dvlbl = [30., 20., 10., 10., 50., 50., 10., 10.]

fig = plt.figure(1, figsize=(10,7))
NPANELS = 8
gs = gridspec.GridSpec(NPANELS/2,2)
axs = []
for i in range(NPANELS):
    axs.append(plt.subplot(gs[i]))

for i in range(len(bs16['ID'])):
    ms = 18
    idx = bs16['ID'][i]
    vinit = pi / 8.0 * ana['vs'][i] # Assumption of vinit    
    if(idx > NPANELS-1): continue
    axs[idx].plot(bs16['t75'][i], bs16['v75'][i], color="grey", marker=symlst[0], markersize=ms)
    axs[idx].plot(bs16['t50'][i], bs16['v50'][i], color="grey", marker=symlst[1], markersize=ms)
    axs[idx].plot(bs16['t25'][i], bs16['v25'][i], color="grey", marker=symlst[2], markersize=ms)

for i in range(len(ana['ID'])):
    if(i > NPANELS-1): continue
    vinit = pi / 8.0 * ana['vs'][i] # Assumption of vinit
    idx = ana['ID'][i]
    if(i == 0 or i == 2): ms = 12
    else: ms = 9
    axs[idx].plot(ana['t75'][i], ana['v75'][i]+vinit, color=clrs[0], marker=symlst[0], markersize=ms)
    axs[idx].plot(ana['t50'][i], ana['v50'][i]+vinit, color=clrs[0], marker=symlst[1], markersize=ms)
    axs[idx].plot(ana['t25'][i], ana['v25'][i]+vinit, color=clrs[0], marker=symlst[2], markersize=ms)
    x = [0.0, ana['t75'][i], ana['t50'][i], ana['t25'][i]]
    y = array([0.0, ana['v75'][i], ana['v50'][i], ana['v25'][i]]) + vinit
    axs[idx].plot(x, y, "-", color=clrs[0])
    # axs[idx].plot([0., 40.], [ana['vs'][i], ana['vs'][i]], "k--")

for i in range(len(ana_sph['ID'])):
    if(i > NPANELS-1): continue
    vinit = pi / 8.0 * ana_sph['vs'][i] # Assumption of vinit
    idx = ana_sph['ID'][i]
    if(i == 0 or i == 2): ms = 12
    else: ms = 9
    axs[idx].plot(ana_sph['t75'][i], ana_sph['v75'][i]+vinit, color=clrs[1], marker=symlst[0], markersize=ms)
    axs[idx].plot(ana_sph['t50'][i], ana_sph['v50'][i]+vinit, color=clrs[1], marker=symlst[1], markersize=ms)
    axs[idx].plot(ana_sph['t25'][i], ana_sph['v25'][i]+vinit, color=clrs[1], marker=symlst[2], markersize=ms)
    x = [0.0, ana_sph['t75'][i], ana_sph['t50'][i], ana_sph['t25'][i]]
    y = array([0.0, ana_sph['v75'][i], ana_sph['v50'][i], ana_sph['v25'][i]]) + vinit
    axs[idx].plot(x, y, "-", color=clrs[1])
    # axs[idx].plot([0., 40.], [ana_sph['vs'][i], ana_sph['vs'][i]], "k--")

if(SHOW_BS16PRED):
    for i in range(len(bs16pred['ID'])):
        idx = bs16pred['ID'][i]    
        axs[idx].plot(bs16pred['tevap'][i], bs16pred['vc'][i], color="black", marker="+", markersize=12)
        print bs16pred['tevap'][i], bs16pred['vc'][i]

for i in range(NPANELS):
    axs[i].text(0.7, 0.1, bs16['Modelname'][i], color='black', fontsize=12, transform=axs[i].transAxes)

axs[0].plot(0.05, 0.9, color="grey", marker=symlst[0], markersize=12, transform=axs[0].transAxes)
axs[0].plot(0.05, 0.75, color="grey", marker=symlst[1], markersize=12, transform=axs[0].transAxes)
axs[0].plot(0.05, 0.60, color="grey", marker=symlst[2], markersize=12, transform=axs[0].transAxes)
axs[0].text(0.08, 0.9, r'$t_{75}$', color='black', fontsize=12, transform=axs[0].transAxes, va='center')
axs[0].text(0.08, 0.75, r'$t_{50}$', color='black', fontsize=12, transform=axs[0].transAxes, va='center')
axs[0].text(0.08, 0.60, r'$t_{25}$', color='black', fontsize=12, transform=axs[0].transAxes, va='center')    

# axs[0].plot(0.80, 0.60, color="blue", marker=symlst[1], markersize=12, transform=axs[0].transAxes)
# axs[0].text(0.83, 0.60, r'$\hat{q}=0.95$', color='blue', fontsize=12, transform=axs[0].transAxes, va='center')
xoff = 0.05
for i in range(NPANELS):
    axs[i].plot(0.80-xoff, 0.9, color="grey", marker=symlst[1], markersize=18, transform=axs[i].transAxes)
    axs[i].text(0.83-xoff, 0.9, "BS16", color='grey', fontsize=12, transform=axs[i].transAxes, va='center')
    axs[i].plot(0.80-xoff, 0.72, color="blue", marker=symlst[1], markersize=12, transform=axs[i].transAxes)   
    axs[i].text(0.83-xoff, 0.72, r'$\hat{q}=0.90$', color='blue', fontsize=12, transform=axs[i].transAxes, va='center')
    axs[i].plot(0.80-xoff, 0.54, color="purple", marker=symlst[1], markersize=12, transform=axs[i].transAxes)       
    axs[i].text(0.83-xoff, 0.54, "spherical", color='purple', fontsize=12, transform=axs[i].transAxes, va='center')

for i in range(NPANELS):
    # axs.append(plt.subplot(gs[i]))
    setp(axs[i].get_xticklabels(),visible=False)
    # setp(axs[-1].get_yticklabels(),visible=False)    
    axs[i].set_xlim(0., 35.)
    axs[i].set_ylim(0., vmax[i])
    axs[i].yaxis.set_ticks(linspace(0., dvlbl[i]*5., 6))
    axs[i].yaxis.set_major_formatter(ticker.FormatStrFormatter("%3d"))
    if((float)(i)/2.0 == (int)(i/2)):
        # setp(axs[i].get_yticklabels(),visible=True)
        axs[i].set_ylabel(r"$\Delta v\ \mathrm{[km/s]}$", fontsize=16)
setp(axs[-1].get_xticklabels(),visible=True)
setp(axs[-2].get_xticklabels(),visible=True)    
setp(axs[-1].get_xticklabels(),visible=True)
axs[-2].xaxis.set_ticks([0, 10, 20, 30])
axs[-1].xaxis.set_ticks([0, 10, 20, 30])    
axs[-2].set_xlabel(r'$t/\tau_{cc}$')
axs[-1].set_xlabel(r'$t/\tau_{cc}$')
fig.subplots_adjust(hspace=0, wspace=0.2, top=0.9, bottom=0.1)

plt.savefig("./figures/vel.pdf")
plt.show()


