from mymod import *
from numpy import genfromtxt, cumsum
from cosmology import tcosmic, acosmic


hparam = 0.7
unit_m = 1.e10 / hparam

CUMULATIVE_DISTRIBUTION = False
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

nbins = 20
aedges = linspace(0., 1., nbins+1)
acen = 0.5 * (aedges[1:] + aedges[:-1])
tcen = (tcosmic(acen) - tcosmic(REDSHIFT)) / 1.e9

acen = 1. / acen - 1. # z
da = 1. / nbins

ZSOLAR = log10(0.0122)

def find_windtime_bin(awind):
    bin_idx = awind / da
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
panels.set_xlabels(r"$z$")
panels.set_ylabels("Fraction")
# panels.set_xlims(0.0, 1.0)
# panels.set_xlims(REDSHIFT, 3.0)
panels.set_xlims(0.0, 10.0)
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
mmin = [11.0, 11.85, 12.85]
mmax = [11.5, 12.15, 13.15]

for modeli in range(len(models)):
    fname = "/home/shuiyao_umass_edu/scidata/gadget3io/"+models[modeli]+"/"+models[modeli]+"_"+zstr+".starinfo"
    soname = "/nas/astro-th-nas/shuiyao/"+models[modeli]+"/so_z"+zstr+".sovcirc"
    foutname = "/home/shuiyao_umass_edu/scidata/gadget3io/"+models[modeli]+"/"+models[modeli]+"_"+zstr+".mprof."
    print ("Doing: ", fname)
    stars = genfromtxt(fname, names=True)
    halos = genfromtxt(soname, names=True)
    halos['Msub'] = log10(halos['Msub'] / 0.7)
    nstars = len(stars)
    stars = stars[stars['WindMass'] > 0.0]
    # stars = stars[stars['a_acc'] > 0.80] # z < 0.25
    stars = stars[stars['a_acc'] > ALIM] # z < 0.25    
    print ("%d out of %d stars selected. " % (len(stars), nstars))
    for mi in range(3):
        mbins_coldw = array([0.0] * nbins)
        mbins_hotw = array([0.0] * nbins)
        for i, s in enumerate(stars):
            hidx = (int)(s['HID'] - 1)
            if(mmin[mi] < halos[hidx]['Msub'] < mmax[mi]):
                aidx = find_windtime_bin(s['WindAge'])
                if(s['Tmax'] < 5.5): mbins_coldw[aidx] += s['WindMass']
                else: mbins_hotw[aidx] += s['WindMass']
        mbins_totw = mbins_coldw + mbins_hotw
        mbins_total = sum(mbins_totw)
        if(CUMULATIVE_DISTRIBUTION):
            mbins_coldw = cumsum(mbins_coldw)
            mbins_hotw = cumsum(mbins_hotw)
            mbins_totw = cumsum(mbins_totw)
        axs[mi].plot(tcen, mbins_coldw/mbins_total, linestyle=lstyles[modeli], color="cyan")
        axs[mi].plot(tcen, mbins_hotw/mbins_total, linestyle=lstyles[modeli], color="magenta")
        axs[mi].plot(tcen, mbins_totw/mbins_total, linestyle=lstyles[modeli], color="orange")
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

plt.savefig("/home/shuiyao_umass_edu/figures/tmp.pdf")
plt.show()

print ("Done.")
