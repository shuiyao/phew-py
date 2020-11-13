import matplotlib.pyplot as plt
from numpy import genfromtxt, log10

models = ["l50n288-phew-m5", "l50n288-phew-m4", "l25n288-phew-m5-spl"]
lstyles = ["--", ":", "-"]
clrs = {'cold':'blue', 'hot':'red', 'igm':'grey'}
Zclrs = {'metalc':'blue', 'metalh':'red', 'metali':'grey'}
fields = ['cold', 'hot', 'igm']
Zfields = ['metalc', 'metalh', 'metali']

from pltastro import frame, draw
import config_mpl
frm = frame.multi(2,1)
pars = frm.params
pars.figsize = (5, 7)
pars.left = 0.2
pars.top = 0.90
pars.bottom = 0.12
panels = frm.panels
panels.set_xlabels(r"$z$")
panels.set_xlims(0.,6.)
panels.ylabels[0] = r"$f_{wind}$"
panels.ylabels[1] = "Z"
fig, axs = draw(frm)

for mi, model in enumerate(models):
    fname = "/home/shuiyao_umass_edu/scidata/" + model + ".wfrac"
    tab = genfromtxt(fname, names=True)
    for field in fields:
        axs[0].plot(tab['z'], tab[field], color = clrs[field], linestyle=lstyles[mi])
    for field in Zfields:
        axs[1].plot(tab['z'], log10(tab[field]), color = Zclrs[field], linestyle=lstyles[mi])
    if(mi == 1): continue
    fname = "/home/shuiyao_umass_edu/scidata/" + model + ".wfracAcc"
    tab = genfromtxt(fname, names=True)
    z, fwind = tab['z'], tab['wind']
    z = 0.5 * (z[0::2] + z[1::2])
    fwind = 0.5 * (fwind[0::2] + fwind[1::2])    
    axs[0].plot(z, fwind, color = "orange", linestyle=lstyles[mi])

from pltastro import legend
lgd = legend.legend(axs[0])
lgd.addLine(("Halo:Cold", "blue", "-", 1))
lgd.addLine(("Halo:Hot", "red", "-", 1))
lgd.addLine(("IGM", "grey", "-", 1))
lgd.addLine(("Accretion", "orange", "-", 1))
lgd.loc = "upper right"
lgd.fontsize = 12
lgd.draw()
lgd = legend.legend(axs[1])
lgd.addLine((models[0], "black", lstyles[0], 1))
lgd.addLine((models[1], "black", lstyles[1], 1))
lgd.addLine((models[2], "black", lstyles[2], 1))
lgd.loc = "upper right"
lgd.fontsize = 12
lgd.draw()

        
plt.savefig("/home/shuiyao_umass_edu/figures/tmp.pdf")
plt.show()
