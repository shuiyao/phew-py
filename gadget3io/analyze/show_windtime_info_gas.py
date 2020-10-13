from mymod import *
from numpy import genfromtxt, cumsum
from cosmology import tcosmic, acosmic


hparam = 0.7
unit_m = 1.e10 / hparam

CUMULATIVE_DISTRIBUTION = False

nbins = 20
aedges = linspace(0., 1., nbins+1)
acen = 0.5 * (aedges[1:] + aedges[:-1])

acen = 1. / acen - 1.
da = 1. / nbins

ZSOLAR = log10(0.0122)

def find_windtime_bin(twind):
    if(twind > 1.34e10): awind = 0.99
    elif(twind < 1.e6): awind = 0.0
    else: awind = acosmic(twind)
    bin_idx = awind / da
    return (int)(bin_idx)

models = ["l25n144-phew-m5-spl", "l25n288-phew-m5-spl"]
lgds = ["25/144,Split", "25/288,Split"]
# models = ["l25n144-phew-m5", "l25n144-phew-m5-spl"]
# lgds = ["PhEW,25/144", "PhEW,25/144,Split"]

lstyles = ["--", "-"]
REDSHIFT = 0.25
if(REDSHIFT == 0.0): zstr = "108"
if(REDSHIFT == 1.0): zstr = "078"
if(REDSHIFT == 2.0): zstr = "058"
if(REDSHIFT == 4.0): zstr = "033"
if(REDSHIFT == 0.25): zstr = "098"
ALIM = acosmic(tcosmic(1./(1.+REDSHIFT)) - 1.e9)

from pltastro import frame, draw
import config_mpl
frm = frame.multi(3,1)
pars = frm.params
pars.figsize = (5, 9)
pars.left = 0.2
pars.top = 0.92
pars.bottom = 0.2
panels = frm.panels
panels.set_xlabels(r"$z$")
panels.set_ylabels("Fraction")
# panels.set_xlims(0.0, 1.0)
panels.set_xlims(REDSHIFT, 3.0)
if(CUMULATIVE_DISTRIBUTION):
    panels.set_ylims(0.0, 1.0)
else:
    panels.set_ylims(0.0, 0.4)
    panels.set_yticks([0.0, 0.1, 0.2, 0.3])

fig, axs = draw(frm)

captions = [
    r"$11.0 < M_\mathrm{vir} < 11.5$",\
    r"$11.85 < M_\mathrm{vir} < 12.15$",\
    r"$12.85 < M_\mathrm{vir} < 13.15$"\    
]
for modeli in range(len(models)):
    for mi, mstr in enumerate(["mh11", "mh12", "mh13"]):
        fname = "/scratch/shuiyao/scidata/gadget3io/"+models[modeli]+"/"+models[modeli]+"_"+zstr+".gas."+mstr
        fauxname = "/scratch/shuiyao/scidata/gadget3io/"+models[modeli]+"/"+models[modeli]+"_"+zstr+".gasaux."+mstr
        print "Doing: ", fname
        gas = genfromtxt(fname, names=True)
        gasaux = genfromtxt(fauxname, names=True)
        gasaux = gasaux[gas['WMass'] > 0.0]        
        gas = gas[gas['WMass'] > 0.0]
        gasaux = gasaux[gas['dr']/gas['Rvir'] < 0.25]
        gas = gas[gas['dr']/gas['Rvir'] < 0.25]
        mbins_coldw = array([0.0] * nbins)
        mbins_hotw = array([0.0] * nbins)
        for i, g in enumerate(gas):
            aidx = find_windtime_bin(gasaux[i]['t'])
            if(g['Tmax'] < 5.5): mbins_coldw[aidx] += g['WMass']
            else: mbins_hotw[aidx] += g['WMass']
        mbins_totw = mbins_coldw + mbins_hotw
        mbins_total = sum(mbins_totw)
        if(CUMULATIVE_DISTRIBUTION):
            mbins_coldw = cumsum(mbins_coldw)
            mbins_hotw = cumsum(mbins_hotw)
            mbins_totw = cumsum(mbins_totw)
        axs[mi].plot(acen, mbins_coldw/mbins_total, linestyle=lstyles[modeli], color="cyan")
        axs[mi].plot(acen, mbins_hotw/mbins_total, linestyle=lstyles[modeli], color="magenta")
        axs[mi].plot(acen, mbins_totw/mbins_total, linestyle=lstyles[modeli], color="orange")
        axs[mi].text(0.05, 0.85, captions[mi], color="black", transform=axs[mi].transAxes)
        
axs[0].set_title("Wind Time, z="+str(REDSHIFT)[:4])

from pltastro import legend
lgd1 = legend.legend(axs[0])
lgd1.loc = "lower right"
lgd1.addLine((lgds[0], "black", lstyles[0], 1))
lgd1.addLine((lgds[1], "black", lstyles[1], 1))
lgd1.draw()
lgd2 = legend.legend(axs[2])
lgd2.loc = "lower right"
lgd2.addLine(("cold wind", "cyan", "-", 1))
lgd2.addLine(("hot wind", "magenta", "-", 1))
lgd2.addLine(("total wind", "orange", "-", 1))
lgd2.draw()

plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
plt.show()

print "Done."
