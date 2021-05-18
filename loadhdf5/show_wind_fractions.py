from myinit import *
import matplotlib.pyplot as plt
from numpy import genfromtxt, log10

models = ["l50n288-phew-m5", "l50n288-phew-m4", "l25n288-phew-m5"]
lboxs = [50., 50., 25.]
lstyles = ["--", ":", "-"]
mclrs = {'mcold':'blue', 'mhot':'red', 'migm':'grey'}
clrs = {'cold':'blue', 'hot':'red', 'igm':'grey'}
Zclrs = {'metalc':'blue', 'metalh':'red', 'metali':'grey'}
fields = ['cold', 'hot', 'igm']
mfields = ['mcold', 'mhot', 'migm']
Zfields = ['metalc', 'metalh', 'metali']

hparam = 0.7
omegab = 0.045
MassInMpc = ac.rhobar * omegab * (ac.mpc/hparam)**3 / ac.msolar / 1.e10 * hparam

from pltastro import frame, draw
import config_mpl
frm = frame.multi(3,1)
pars = frm.params
pars.figsize = (5, 9)
pars.left = 0.2
pars.top = 0.90
pars.bottom = 0.12
panels = frm.panels
panels.set_xlabels(r"$z$")
panels.set_xlims(0.,10.)
panels.ylabels[0] = r"$f_{phase}$"
panels.ylabels[1] = r"$f_{wind}$"
panels.ylabels[2] = "Z"
panels.ylims[0] = (0.0, 1.0)
panels.ylims[1] = (0.0, 0.6)
panels.ylims[2] = (-5., -2.)
panels.yticks[0] = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
panels.yticks[1] = [0.0, 0.2, 0.4, 0.6]
panels.yticks[2] = [-5., -4., -3.]
fig, axs = draw(frm)

for mi, model in enumerate(models):
    fname = "/home/shuiyao_umass_edu/scidata/" + model + "/" + model + ".wfrac"
    tab = genfromtxt(fname, names=True)
    for field in fields:
        axs[1].plot(tab['z'], tab[field], color = clrs[field], linestyle=lstyles[mi])
    for field in Zfields:
        axs[2].plot(tab['z'], log10(tab[field]), color = Zclrs[field], linestyle=lstyles[mi])
    MassTotal = MassInMpc * (lboxs[mi]) ** 3
    for field in mfields:
        axs[0].plot(tab['z'], tab[field]/MassTotal, color = mclrs[field], linestyle=lstyles[mi])
    fname = "/home/shuiyao_umass_edu/scidata/" + model + "/" + model + ".wfracAcc"
    tab = genfromtxt(fname, names=True)
    z, fwind, metal = tab['z'], tab['wind'], tab['metals']
    # z = 0.5 * (z[0::2] + z[1::2])
    # fwind = 0.5 * (fwind[0::2] + fwind[1::2])
    metal = log10(metal)
    # metal = 0.5 * (metal[0::2] + metal[1::2])
    axs[1].plot(z, fwind, color = "orange", linestyle=lstyles[mi])
    axs[2].plot(z, metal, color = "orange", linestyle=lstyles[mi])    
    

from pltastro import legend
lgd = legend.legend(axs[1])
lgd.addLine(("Halo:Cold", "blue", "-", 1))
lgd.addLine(("Halo:Hot", "red", "-", 1))
lgd.addLine(("IGM", "grey", "-", 1))
lgd.addLine(("Accretion", "orange", "-", 1))
lgd.loc = "upper right"
lgd.fontsize = 12
lgd.draw()
lgd = legend.legend(axs[0])
lgd.addLine((models[0], "black", lstyles[0], 1))
lgd.addLine((models[1], "black", lstyles[1], 1))
lgd.addLine((models[2], "black", lstyles[2], 1))
lgd.loc = "upper right"
lgd.fontsize = 12
lgd.draw()

        
plt.savefig("/home/shuiyao_umass_edu/figures/tmp.pdf")
plt.show()
