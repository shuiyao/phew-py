from myinit import *
from numpy import genfromtxt, cumsum
from cosmology import tcosmic, acosmic
import matplotlib.pyplot as plt

hparam = 0.7
unit_m = 1.e10 / hparam

nbins = 40
SIGMAX = 150.0
sigedges = linspace(0., SIGMAX, nbins+1)
sigcen = 0.5 * (sigedges[1:] + sigedges[:-1])
dsig = SIGMAX / nbins

ZSOLAR = log10(0.0122)

def find_siggal_bin(sig):
    if(sig >= SIGMAX): sig = SIGMAX - 1.0
    bin_idx = sig / dsig
    return (int)(bin_idx)

# models = ["l50n288-phew-m5", "l50n576-phew-m5"]
# lgds = ["l50n288-phew-m5", "l50n576-phew-m5"]
models = ["l25n144-phew-m5-spl", "l25n288-phew-m5"]
lgds = ["l25n144-phew-m5-spl", "l25n288-phew-m5"]

lstyles = ["--", "-"]
REDSHIFT = 0.0
if(REDSHIFT == 0.0): zstr = "108"
if(REDSHIFT == 1.0): zstr = "078"
if(REDSHIFT == 2.0): zstr = "058"
if(REDSHIFT == 4.0): zstr = "033"
if(REDSHIFT == 0.25): zstr = "098"
ALIM = acosmic(tcosmic(1./(1.+REDSHIFT)) - 1.e9)
ALIm = 0

from pltastro import frame, draw
import config_mpl
frm = frame.multi(2,1)
pars = frm.params
pars.figsize = (5, 7)
pars.left = 0.2
pars.top = 0.92
pars.bottom = 0.2
panels = frm.panels
panels.set_xlabels(r"$\sigma_{gal}$")
panels.set_xlims(0.0, SIGMAX)
panels.set_ylabels("Mass")
# panels.set_ylabels("Fraction")
# panels.set_ylims(0.0, 1.0)

fig, axs = draw(frm)

captions = [
    r"$11.0 < M_\mathrm{vir} < 11.5$",\
    r"$11.85 < M_\mathrm{vir} < 12.15$",\
    r"$12.85 < M_\mathrm{vir} < 13.15$"\    
]
mmin = [11.0, 11.85, 12.85]
mmax = [11.5, 12.15, 13.15]

for modeli in range(len(models)):
    fname = DIRS['SCIDATA']+models[modeli]+"/"+models[modeli]+"_"+zstr+".starinfo"
    soname = DIRS['DATA']+models[modeli]+"/so_z"+zstr+".sovcirc"
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
        mbins_siggal = array([0.0] * nbins)
        for i, s in enumerate(stars):
            hidx = (int)(s['HID'] - 1)
            if(mmin[mi] < halos[hidx]['Msub'] < mmax[mi]):
                aidx = find_siggal_bin(s['WindSig'])
                mbins_siggal[aidx] += s['WindMass']
        # mbins_siggal = cumsum(mbins_siggal)
        # axs[mi].plot(sigcen, mbins_siggal/mbins_siggal[-1], linestyle=lstyles[modeli], color="black")
        axs[mi].plot(sigcen, mbins_siggal, linestyle=lstyles[modeli], color="black")        
        axs[mi].text(0.55, 0.85, captions[mi], color="black", transform=axs[mi].transAxes)
        
#axs[0].set_title("$\sigma_{gal}$, z="+str(REDSHIFT)[:4])

from pltastro import legend
lgd1 = legend.legend(axs[0])
lgd1.loc = "lower right"
lgd1.size = 8
lgd1.addLine((lgds[0], "black", lstyles[0], 1))
lgd1.addLine((lgds[1], "black", lstyles[1], 1))
lgd1.draw()

plt.savefig(DIRS['FIGURE']+"tmp.pdf")
plt.show()

print ("Done.")
