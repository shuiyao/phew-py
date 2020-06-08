from mymod import *
from numpy import genfromtxt

fname = "/scratch/shuiyao/scidata/gadget3io/l25n144-phew-fa/l25n144-phew-fa_078.gas.mh13"

hparam = 0.7
unit_m = 1.e10 / hparam

nbins = 30
redges = linspace(0., 1., nbins+1)
dr = 1. / nbins

def find_radial_bin(r):
    if(r > 1.0): return -1
    bin_idx = r / dr
    return (int)(bin_idx)

from pltastro import frame, draw
import config_mpl
frm = frame.multi(3,1)
pars = frm.params
pars.figsize = (5, 9)
pars.left = 0.2
pars.top = 0.92
pars.bottom = 0.2
panels = frm.panels
panels.set_xlabels(r"$R/R_\mathrm{vir}$")
panels.set_ylabels("")
# panels.ylabels[1] = "Mass [arbitrary units]"
panels.ylabels[1] = "Mass Fraction"

fig, axs = draw(frm)

captions = [
    r"$11.0 < M_\mathrm{vir} < 11.5$",\
    r"$11.85 < M_\mathrm{vir} < 12.15$",\
    r"$12.85 < M_\mathrm{vir} < 13.15$"\    
]

for mi, mstr in enumerate(["mh11", "mh12", "mh13"]):
    fname = "/scratch/shuiyao/scidata/gadget3io/l25n144-phew-fa/l25n144-phew-fa_078.gas." + mstr
    # fname = "/scratch/shuiyao/scidata/gadget3io/l25n144-phewoff/l25n144-phewoff_078.gas." + mstr    

    tab = genfromtxt(fname, names=True)

    mass = array([0.0] * (nbins + 1))
    mwind = array([0.0] * (nbins + 1))
    mmix = array([0.0] * (nbins + 1))

    for i in range(len(tab)):
        part = tab[i]
        bidx = find_radial_bin(part['dr'] / part['Rvir'])
        mass[bidx] += part['Mass']
        mmix[bidx] += part['WMass']
        if(part['Mc'] > 0):
            mwind[bidx] += part['Mass']
            mmix[bidx] -= part['WMass']
        if(part['Mc'] < 0):            
            mwind[bidx] += part['Mass'] - part['WMass']
    fwind = mwind/mass

    axs[mi].plot(redges, mass/mass, "k-")
    axs[mi].plot(redges, mwind/mass, "g-")
    axs[mi].plot(redges, mmix/mass, "-", color="cyan")
    axs[mi].text(0.2, 0.8, captions[mi], transform=axs[mi].transAxes, fontsize=16)
axs[0].set_title("rec-fast, z=1")


for mi, mstr in enumerate(["mh11", "mh12", "mh13"]):
    fname = "/scratch/shuiyao/scidata/gadget3io/l25n144-phewoff/l25n144-phewoff_078.gas." + mstr    

    tab = genfromtxt(fname, names=True)

    mass = array([0.0] * (nbins + 1))
    mwind = array([0.0] * (nbins + 1))
    mmix = array([0.0] * (nbins + 1))

    for i in range(len(tab)):
        part = tab[i]
        bidx = find_radial_bin(part['dr'] / part['Rvir'])
        mass[bidx] += part['Mass']
        if(part['Mc'] < 0):            
            mwind[bidx] += part['Mass']
    fwind = mwind/mass

    axs[mi].plot(redges, mass/mass, "k--")
    axs[mi].plot(redges, mwind/mass, "g--")
    # axs[mi].text(0.2, 0.8, captions[mi], transform=axs[mi].transAxes, fontsize=16)

# plt.plot(redges, mass, "k-")
# plt.plot(redges, mwind, "g-")
# plt.plot(redges, mmix, "-", color="cyan")
plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
plt.show()

print "Done."
