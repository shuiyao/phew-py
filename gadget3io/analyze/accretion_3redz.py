from myinit import *
import pandas as pd
from cosmology import tcosmic, acosmic
import matplotlib.pyplot as plt

# Input: .starinfo

models = ["l50n288-phewoff", "l50n288-phew-m5", "l50n576-phew-m5"]
lstyles = {"l50n288-phewoff":":", "l50n288-phew-m5":"--", "l50n576-phew-m5":"-"}
#model = "l50n576-phew-m5"
cols_halos = ["HID", "Mvir", "Rvir", "R_vmax", "Vmax", "Npart", "Msub", "Rsub"] 
dtypes_halos = {'HID':'int', 'Npart':'int'}  
FEATURES = ['a_acc', 'Mass', 'WindMass', 'Tmax', 'Z', 'HID']
REDSHIFT = 1.0
HPARAM = 0.7
TCUT = 5.5

def find_ascale_range_for_redshift(z, dtime=1.e9):
    amax = 1./(z+1.)
    amin = acosmic(tcosmic(amax) - dtime)
    return amin, amax

def draw_for_one_redshift(tab, halos, zred, ax=None, lstyle="-"):
    if(zred == -1): amin, amax = 0.0, 1.0
    else:
        amin, amax = find_ascale_range_for_redshift(zred)
    select = str(amin)+' < a_acc < '+ str(amax) + ' & GID == HID & a_last <= 0.0'
    tabz1 = tab.query(select)[FEATURES]
    tabz1 = tabz1.join(halos[['HID','Msub']].set_index('HID'), on=tabz1['HID'])
    tabz1['Msub'] = log10(tabz1['Msub'] / HPARAM)
    mbins = linspace(11., 14.0, 21)
    mcens = 0.5 * (mbins[1:] + mbins[:-1])
    tabz1['Mbin'] = pd.cut(tabz1['Msub'], bins=mbins, labels=mcens)
    groups = tabz1.groupby(['Mbin', tabz1['Tmax'] < TCUT]) # True: cold
    mtab = groups.sum().unstack()
    logmsub = mtab.index.values
    mhot = mtab['Mass', False]
    mcold = mtab['Mass', True]
    mwhot = mtab['WindMass', False]
    mwcold = mtab['WindMass', True]
    mtot = mhot + mcold
    mhot = mhot - mwhot
    mcold = mcold - mwcold
    mwind = mwhot + mwcold
    macc = {'cold':mcold, 'hot':mhot, 'wind':mwind, 'wcold':mwcold, 'whot':mwhot, 'total':mtot}

    ax.plot(logmsub, mcold/mtot, "b", linestyle=lstyle)
    ax.plot(logmsub, mhot/mtot, "r", linestyle=lstyle)
    ax.plot(logmsub, mwind/mtot, "g", linestyle=lstyle)

fig, axs = plt.subplots(3, 1, figsize=(5,10))
for model in models:
    print("Plotting for model: ", model)
    data_sample = DIRS['SCIDATA'] + model + "/" + model + "_108.starinfo"
    halo_sample = DIRS['DATA'] + model + "/" + "so_z108.sovcirc"

    halos = pd.read_csv(halo_sample, sep='\s+', skiprows=1, names=cols_halos, dtype=dtypes_halos)
    tab = pd.read_csv(data_sample, sep='\s+', header=0)
    # zred = REDSHIFT
    draw_for_one_redshift(tab=tab, halos=halos, zred=-1, ax=axs[0], lstyle=lstyles[model])
    draw_for_one_redshift(tab=tab, halos=halos, zred=0.0, ax=axs[1], lstyle=lstyles[model])
    draw_for_one_redshift(tab=tab, halos=halos, zred=2.0, ax=axs[2], lstyle=lstyles[model])

redshifts = ["All", "0.0", "2.0"]
for i, ax in enumerate(axs):
    ax.set_xlim(11.0, 14.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xticks([])
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])    
    ax.set_title("Accretion at z ~ %s" % (redshifts[i]))
    ax.set_ylabel("Fraction")
axs[2].set_xticks([11.0, 12.0, 13.0, 14.0])
axs[2].set_xlabel(r"$\log(M_{vir}/M_\odot)$")

from pltastro import legend
lgd = legend.legend(ax=axs[0])
lgd.addLine(('cold', 'blue', '-', 1))
lgd.addLine(('hot', 'red', '-', 1))
lgd.addLine(('wind', 'green', '-', 1))
lgd.loc = 'upper right'
lgd.draw()
lgd2 = legend.legend(ax=axs[2])
for model in models: lgd2.addLine((model, 'black', lstyles[model], 1))
lgd2.loc = 'upper right'
lgd2.draw()

plt.savefig(DIRS['FIGURE']+"tmp.pdf")
plt.show()

print ("Done.")
