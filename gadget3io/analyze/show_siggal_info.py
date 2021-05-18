from myinit import *
from numpy import genfromtxt, cumsum
from cosmology import tcosmic, acosmic
import matplotlib.pyplot as plt
import pandas as pd

hparam = 0.7
unit_m = 1.e10 / hparam

nbins = 40
SIGMAX = 150.0
sigbins = linspace(0., SIGMAX, nbins+1)
sigcen = 0.5 * (sigbins[1:] + sigbins[:-1])
dsig = SIGMAX / nbins

ZSOLAR = log10(0.0122)

def find_siggal_bin(sig):
    if(sig >= SIGMAX): sig = SIGMAX - 1.0
    bin_idx = sig / dsig
    return (int)(bin_idx)

models = ["l50n288-phew-m5", "l50n576-phew-m5"]
lgds = ["l50n288-phew-m5", "l50n576-phew-m5"]
ncpus = [256, 1024]
# models = ["l25n144-phew-m5-spl", "l25n288-phew-m5"]
# lgds = ["l25n144-phew-m5-spl", "l25n288-phew-m5"]
# ncpus = [128, 256]

lstyles = ["--", "-"]
REDSHIFT = 0.0
if(REDSHIFT == 0.0): zstr = "108"
if(REDSHIFT == 1.0): zstr = "078"
if(REDSHIFT == 2.0): zstr = "058"
if(REDSHIFT == 4.0): zstr = "033"
if(REDSHIFT == 0.25): zstr = "098"
ALIM = acosmic(tcosmic(1./(1.+REDSHIFT)) - 1.e9)
# ALIM = 0

from pltastro import frame, draw
import config_mpl
frm = frame.multi(2,1)
pars = frm.params
pars.figsize = (5, 7)
pars.left = 0.2
pars.top = 0.90
pars.bottom = 0.1
panels = frm.panels
panels.set_xlabels(r"$\sigma_{gal} [km/s]$")
panels.set_xlims(0.0, SIGMAX)
panels.set_ylabels("Mass") # <--------------------------------
# panels.set_ylabels("Fraction")
# panels.set_ylims(0.0, 1.0)

fig, axs = draw(frm)

LBOX_MPC = 50.0
# Upper Panel:
zstrs = ["z2", "z0"]
clrs = ["blue", "red"]
dtatz = [tcosmic(0.35) - tcosmic(0.333333),\
         tcosmic(0.85) - tcosmic(0.833333)]

cols = ["atime", "ID", "mass", "siggal", "mgal", "v", "mstar", "mgas", "sfr", 'zcool', 'x', 'y', 'z']

def draw_init_siggal(fname, ax, dt, lstyle="-", clr="black"):
    df = pd.read_csv(fname, sep='\s+', names=cols)
    df['sigbin'] = pd.cut(df['siggal'], bins=sigbins, labels=sigcen)
    grps = df.groupby('sigbin')
    mtot = grps['mass'].sum()
    mtot = mtot * 1.e10/0.7/(LBOX_MPC * dt)
    ax.plot(sigcen, mtot, linestyle=lstyle, color=clr)

print("Drawing for ", models[0])    
for zi, zstr in enumerate(zstrs):
    fname = DIRS['DATA']+models[0]+"/WINDS/winds."+zstr
    draw_init_siggal(fname, axs[0], dtatz[zi], "--", clrs[zi])
print("Drawing for ", models[1])        
for zi, zstr in enumerate(zstrs):
    fname = DIRS['DATA']+models[1]+"/WINDS/winds."+zstr
    draw_init_siggal(fname, axs[0], dtatz[zi], "-", clrs[zi])

axs[0].set_ylabel(r"$\dot{M}_{w} [M_\odot (cMpc)^{-1} yr^{-1}]$")
# axs[0].set_ylims(0.0, 2.5)
# axs[0].set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])

# Bottom Panel: accretion
zstr = "108"
for modeli in range(len(models)):
    outname = DIRS['SCIDATA']+models[modeli]+"/"+"siggal_accretion_z0.csv"
    if(os.path.exists(outname)):
        print("Reading {}".format(outname))
        df = pd.read_csv(outname)
    else:
        fname = DIRS['SCIDATA']+models[modeli]+"/"+models[modeli]+"_"+zstr+".starinfo"
        soname = DIRS['DATA']+models[modeli]+"/so_z"+zstr+".sovcirc"
        print ("Doing: ", fname)
        stars = genfromtxt(fname, names=True)
        # halos = genfromtxt(soname, names=True)
        # halos['Msub'] = log10(halos['Msub'] / 0.7)

        nstars = len(stars)
        stars = stars[stars['WindMass'] > 0.0]
        # stars = stars[stars['a_acc'] > 0.80] # z < 0.25
        stars = stars[stars['a_acc'] > ALIM] # z < 0.25    
        print ("%d out of %d stars selected. " % (len(stars), nstars))
        mbins_siggal = array([0.0] * nbins)
        for i, s in enumerate(stars):
            hidx = (int)(s['HID'] - 1)
            aidx = find_siggal_bin(s['WindSig'])
            mbins_siggal[aidx] += s['WindMass']
        mbins_siggal *= 1./(1.989e33 * LBOX_MPC * 1.e9)
        df = pd.DataFrame({'sigcen':sigcen, 'mass':mbins_siggal})
        print("Writing {}".format(outname))    
        df.to_csv(outname)
    
    axs[1].plot(df.sigcen, df.mass, linestyle=lstyles[modeli], color="black")

axs[1].set_ylabel(r"$\dot{M}_{acc} [M_\odot (cMpc)^{-1} yr^{-1}]$")
axs[1].set_ylim(0.0, 2.5)
axs[1].set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])

from pltastro import legend
lgd1 = legend.legend(axs[0])
lgd1.loc = "upper right"
lgd1.size = 8
lgd1.addLine(("z=2", "blue", "-", 1))
lgd1.addLine(("z=0", "red", "-", 1))
lgd1.draw()

lgd1 = legend.legend(axs[1])
lgd1.loc = "upper right"
lgd1.size = 8
lgd1.addLine((lgds[0], "black", lstyles[0], 1))
lgd1.addLine((lgds[1], "black", lstyles[1], 1))
lgd1.draw()

plt.savefig(DIRS['FIGURE']+"tmp.pdf")
plt.show()

print ("Done.")
