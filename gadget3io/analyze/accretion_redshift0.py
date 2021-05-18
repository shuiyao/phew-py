from myinit import *
import pandas as pd
from cosmology import tcosmic, acosmic
import matplotlib.pyplot as plt

# Input: .starinfo

# models = ["l50n288-phewoff", "l50n288-phew-m5", "l50n576-phew-m5"]
models = ["l50n288-phewoff", "l50n576-phew-m5"]
lstyles = {"l50n288-phewoff":":", "l50n288-phew-m5":"--", "l50n576-phew-m5":"-"}
#model = "l50n576-phew-m5"
cols_halos = ["HID", "Mvir", "Rvir", "R_vmax", "Vmax", "Npart", "Msub", "Rsub"] 
dtypes_halos = {'HID':'int', 'Npart':'int'}  
FEATURES = ['a_acc', 'Mass', 'StarMass', 'WindMass', 'Tmax', 'Z', 'HID']
REDSHIFT = 1.0
HPARAM = 0.7
TCUT = 5.5

def find_ascale_range_for_redshift(z, dtime=1.e9):
    amax = 1./(z+1.)
    amin = acosmic(tcosmic(amax) - dtime)
    return amin, amax

def draw_for_one_redshift(tab, halos, axs, lstyle="-"):
    amin, amax = 0.0, 1.0
    select = str(amin)+' < a_acc < '+ str(amax) + ' & GID == HID & a_last <= 0.0'
    tabz1 = tab.query(select)[FEATURES]
    tabz1 = tabz1.join(halos[['HID','Msub']].set_index('HID'), on=tabz1['HID'])
    tabz1['LogMsub'] = log10(tabz1['Msub'] / HPARAM)
    tabz1['WindMass'] *= (tabz1['StarMass'] / tabz1['Mass'])
    mbins = linspace(11., 14.0, 21)
    mcens = 0.5 * (mbins[1:] + mbins[:-1])
    tabz1['Mbin'] = pd.cut(tabz1['LogMsub'], bins=mbins, labels=mcens)
    groups = tabz1.groupby(['Mbin', tabz1['Tmax'] < TCUT]) # True: cold
    mtotsub = tabz1.drop_duplicates('Msub').groupby('Mbin')['Msub'].sum()
    mtotsub *= 2.e33 * (0.045 / 0.3) / HPARAM
    mtab = groups.sum().unstack()
    logmsub = mtab.index.values
    mhot = mtab['StarMass', False]
    mcold = mtab['StarMass', True]
    mwhot = mtab['WindMass', False]
    mwcold = mtab['WindMass', True]
    mtot = mhot + mcold
    mhot = mhot - mwhot
    mcold = mcold - mwcold
    mwind = mwhot + mwcold
    macc = {'cold':mcold, 'hot':mhot, 'wind':mwind, 'wcold':mwcold, 'whot':mwhot, 'total':mtot}

    # FBEHROOZI = DIRS['SCI']+"REFERENCES/behroozi13/smmr/"
    # snapnum, zstr_behroozi = "108", "0.10"
    # mh, ms, err1, err2 = ioformat.rcol(FBEHROOZI+"c_smmr_z"+zstr_behroozi+"_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat", [0,1,2,3], linestart=1)
    # ms = array(ms) - log10(0.156)
    # axs[0].plot(mh, ms, "o", color="purple", markersize=6)
    # axs[0].errorbar(mh, ms, yerr=[err1, err2], color="purple")

    axs[0].plot(logmsub, log10(mcold/mtotsub), "b", linestyle=lstyle, linewidth=2)
    axs[0].plot(logmsub, log10(mhot/mtotsub), "r", linestyle=lstyle, linewidth=2)
    axs[0].plot(logmsub, log10(mwind/mtotsub), "g", linestyle=lstyle, linewidth=2)
    axs[0].plot(logmsub, log10(mtot/mtotsub), "k", linestyle=lstyle, linewidth=2)    

    axs[1].plot(logmsub, mcold/mtot, "b", linestyle=lstyle, linewidth=2)
    axs[1].plot(logmsub, mhot/mtot, "r", linestyle=lstyle, linewidth=2)
    axs[1].plot(logmsub, mwind/mtot, "g", linestyle=lstyle, linewidth=2)

    # axs[0].plot([10.9, 10.9], [-4.0, -0.5], "k--")
    # axs[0].plot([10.0, 10.0], [-4.0, -0.5], "k-")
    # axs[1].plot([10.9, 10.9], [0.0, 1.1], "k--")
    # axs[1].plot([10.0, 10.0], [0.0, 1.1], "k-")

from pltastro import *
import pltastro
frm = pltastro.frame.multi(2, 1)
frm.panels.ylabels[0] = r"$\log(\frac{M_*}{(\Omega_b/\Omega_m)M_{h}})$"
frm.panels.ylabels[1] = r"$f_{acc}$"
frm.panels.set_xlabels(r"$\log(M_{vir}/M_\odot)$")
frm.panels.set_xlims(10., 14.)
frm.panels.ylims[1] = (0.0, 1.1)
frm.panels.yticks[1] = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
frm.panels.ylims[0] = (-3.0, -0.0)
frm.panels.yticks[0] = [-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0]
frm.panels.ytickformat = ['%3.1f', "%3.1f"]
frm.panels.set_xlims(11.0, 13.5)
frm.params.figsize = (6, 7)
frm.params.height_ratios = [2, 1]
frm.params.left = 0.20
fig, axs = draw(frm)

for model in models:
    print("Plotting for model: ", model)
    data_sample = DIRS['SCIDATA'] + model + "/" + model + "_108.starinfo"
    halo_sample = DIRS['DATA'] + model + "/" + "so_z108.sovcirc"

    halos = pd.read_csv(halo_sample, sep='\s+', skiprows=1, names=cols_halos, dtype=dtypes_halos)
    tab = pd.read_csv(data_sample, sep='\s+', header=0)
    # zred = REDSHIFT
    draw_for_one_redshift(tab=tab, halos=halos, axs=axs, lstyle=lstyles[model])

from pltastro import legend
lgd = legend.legend(ax=axs[1])
lgd.addLine(('cold', 'blue', '-', 1))
lgd.addLine(('hot', 'red', '-', 1))
lgd.addLine(('wind', 'green', '-', 1))
lgd.loc = 'upper right'
lgd.draw()
lgd2 = legend.legend(ax=axs[0])
for model in models: lgd2.addLine((model, 'black', lstyles[model], 1))
lgd2.loc = 'lower right'
lgd2.draw()

plt.savefig(DIRS['FIGURE']+"tmp.pdf")
plt.show()

print ("Done.")
