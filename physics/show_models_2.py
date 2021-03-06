import ioformat
import matplotlib.pyplot as plt
from scipy import sqrt, array, log10, linspace
from astropy.table import Table
from astropy.io import ascii
import cosmology
from matplotlib import gridspec
from pylab import setp

bs16 = Table.read("models_bs16.dat", format="ascii", data_start=0)
bs16pred = Table.read("models_bs16_predictions.dat", format="ascii", data_start=0)
# ana = Table.read("models_modela.dat", format="ascii", data_start=0)
ana = Table.read("models_180727.dat", format="ascii", data_start=0)

symlst = ["*", "^", ".", ".", ".", "."]
clrs = ["blue", "purple", "red", "green", "yellow", "purple", "brown"]
vmax = [120., 120., 250., 80., 210., 320.]
dvlbl = [20., 20., 50., 15., 40., 60.]

fig = plt.figure(1, figsize=(6,9))
gs = gridspec.GridSpec(3,1,height_ratios=[1]*3)
axs = []
for i in range(3):
    axs.append(plt.subplot(gs[i]))
    setp(axs[-1].get_xticklabels(),visible=False)
    axs[-1].set_xlim(0., 40.)
    axs[-1].set_ylim(0., vmax[i])
    axs[-1].yaxis.set_ticks(linspace(0., dvlbl[i]*5., 6))
    axs[-1].set_ylabel("dv")
setp(axs[-1].get_xticklabels(),visible=True)
axs[-1].set_xlabel("t/tcc")
fig.subplots_adjust(hspace=0, wspace=0, top=0.9, bottom=0.1)

for i in range(len(ana['ID'])):
    idx = ana['ID'][i]
    if(i == 0 or i == 3): ms = 18
    else: ms = 9
    axs[0].plot(ana['t75'][i], ana['v75'][i], color=clrs[ana['Group'][i]], marker=symlst[0], markersize=ms)
    axs[1].plot(ana['t50'][i], ana['v50'][i], color=clrs[ana['Group'][i]], marker=symlst[1], markersize=ms)
    axs[2].plot(ana['t25'][i], ana['v25'][i], color=clrs[ana['Group'][i]], marker=symlst[2], markersize=ms)

for i in range(len(bs16['ID'])):
    ms = 12
    idx = bs16['ID'][i]
    if(idx > 5): continue
    axs[0].plot(bs16['t75'][i], bs16['v75'][i], color="black", marker=symlst[0], markersize=ms)
    axs[1].plot(bs16['t50'][i], bs16['v50'][i], color="black", marker=symlst[1], markersize=ms)
    axs[2].plot(bs16['t25'][i], bs16['v25'][i], color="black", marker=symlst[2], markersize=ms)
    # axs[idx].plot([0., 40.], [bs16['vexp'][i], bs16['vexp'][i]], "k--")

for i in range(len(bs16pred['ID'])):
    idx = bs16pred['ID'][i]
    for j in range(3):
        axs[j].plot(bs16pred['tevap'][i], bs16pred['vc'][i], color="black", marker="+", markersize=12)
    print bs16pred['tevap'][i], bs16pred['vc'][i]

# for i in range(6):
#     axs[i].text(0.7, 0.1, bs16['Modelname'][i], color='black', fontsize=12, transform=axs[i].transAxes)
    
plt.show()


