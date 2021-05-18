# Input: gasprof_*

from myinit import *
import matplotlib.pyplot as plt

hparam = 0.7
unit_m = 1.e10 / hparam

nbins = 30
redges = linspace(0., 1., nbins+1)
redges = redges[:-1]
dr = 1. / nbins

ZSOLAR = log10(0.0122)

#MODE = "Metal"
MODE = "Mass"

def find_radial_bin(r):
    if(r > 1.0): return -1
    bin_idx = r / dr
    return (int)(bin_idx)

# ---- GIZMO Vs. PhEW
# models = ["l50n288-phewoff", "l50n288-phew-m5", "l50n576-phew-m5"]
# lgds = ["l50n288-phewoff", "l50n288-phew-m5", "l50n576-phew-m5"]
# lstyles = [":", "--", "-"]
models = ["l50n288-phewoff", "l50n576-phew-m5"]
lgds = ["l50n288-phewoff", "l50n576-phew-m5"]
lstyles = [":", "-"]
REDSHIFT = 0.25
# REDSHIFT = 1.0
if(REDSHIFT == 0.0): zstr = "108"
if(REDSHIFT == 0.25): zstr = "098"
if(REDSHIFT == 1.0): zstr = "078"
if(REDSHIFT == 2.0): zstr = "058"
if(REDSHIFT == 4.0): zstr = "033"

from pltastro import frame, draw
import config_mpl
frm = frame.multi(3,2)
pars = frm.params
pars.figsize = (10, 9)
pars.left = 0.1
pars.right = 0.9
pars.top = 0.92
pars.wspace = 0.3
pars.bottom = 0.15
panels = frm.panels
panels.set_xlabels(r"$R/R_\mathrm{vir}$")
panels.set_ylabels("")
# panels.ylabels[1] = "Mass [arbitrary units]"
panels.ylabels[2] = "Metal Fraction"
panels.ylabels[3] = r"$\log(Z/Z_\odot)$"
panels.yticksON = [True] * 6

fig, axs = draw(frm)

captions = [
    r"$11.0 < M_\mathrm{vir} < 11.5$",\
    r"$11.85 < M_\mathrm{vir} < 12.15$",\
    r"$12.85 < M_\mathrm{vir} < 13.15$"\    
]

for modeli in range(len(models)):
    for mi, mstr in enumerate(["mh11", "mh12", "mh13"]):
        foutname = DIRS['SCIDATA'] + models[modeli] + "/gasprof_"+zstr+"_"+mstr
        print ("Plotting: "+ foutname)

        tab = genfromtxt(foutname, names=True)
        tab = tab[:-1]
        mass = tab['Mass']
        mwind = tab['Mwind']
        mwcold = tab['Mwcold']
        mwhot = tab['Mwhot']
        mcold = tab['Mcold']
        mhot = tab['Mhot']
        mism = tab['Mism']
        mmaxhot = tab['Mmaxhot']
        mmaxcold = tab['Mmaxcold']
        mmaxwhot = tab['Mmaxwhot']
        mmaxwcold = tab['Mmaxwcold']
        mzcold = tab['Mzcold']
        mzhot = tab['Mzhot']
        mzism = tab['Mzism']

        if(models[modeli] == "l50n288-phewoff"):
            mzwind = tab['Mzwind']
            ZWIND = log10(mzwind / (mwcold + mwhot))
            mz = mzcold + mzhot + mzism
            mzwind = 0.0
        else:
            mzwind = tab['Mzwind']
            ZWIND = log10(mzwind / mwind)
            mz = mzcold + mzhot + mzwind + mzism
            # mzwind is Z in PhEW particles

        # logT based
        # axs[2*mi].plot(redges, mass/mass, color="black", linestyle=lstyles[modeli])
        axs[2*mi].plot(redges, mzcold/mz, color="blue", linestyle=lstyles[modeli])
        axs[2*mi].plot(redges, mzhot/mz, color="red", linestyle=lstyles[modeli])
        axs[2*mi].plot(redges, mzwind/mz, color="green", linestyle=lstyles[modeli])
        # axs[2*mi].plot(redges, mzism/mz, color="brown", linestyle=lstyles[modeli])
        axs[2*mi].set_ylim(0.0, 1.0)
        axs[2*mi].set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])

        # Z
        # PhEW: mzcold and mzhot are mZ in cold and hot gas
        #       mzwind is mZ in active PhEW particles
        # PhEWOFF: 
        
        axs[2*mi+1].plot(redges, log10((mzcold+mzhot)/(mcold+mwcold+mhot+mwhot))-ZSOLAR, color="black", linestyle=lstyles[modeli])        
        axs[2*mi+1].plot(redges, log10(mzcold/(mcold+mwcold))-ZSOLAR, color="blue", linestyle=lstyles[modeli])
        axs[2*mi+1].plot(redges, log10(mzhot/(mhot+mwhot))-ZSOLAR, color="red", linestyle=lstyles[modeli])
        # axs[2*mi+1].plot(redges, log10(mzism/mism)-ZSOLAR, color="brown", linestyle=lstyles[modeli])
        # axs[2*mi+1].plot(redges, ZWIND-ZSOLAR, color="green", linestyle=lstyles[modeli])
        axs[2*mi+1].set_ylim(-1.5, 0.5)
        axs[2*mi+1].set_yticks([-1.5, -1.0, -0.5, 0.0])

        # Text
        axs[2*mi+1].text(0.05, 0.85, captions[mi], transform=axs[2*mi+1].transAxes, fontsize=12)

            
    axs[0].set_title("Metal Profiles, z="+str(REDSHIFT)[:3])

    # axs[1].set_title("Wind Profiles, z="+str(REDSHIFT)[:3])
    axs[1].set_title("Metal Profiles, z="+str(REDSHIFT)[:3])    

    from pltastro import legend
    lgd1 = legend.legend(axs[0])
    lgd1.loc = "upper right"
    for i in range(len(lgds)):
        lgd1.addLine((lgds[i], "black", lstyles[i], 1))
    lgd1.draw()
    lgd2 = legend.legend(axs[4])
    lgd2.loc = "upper right"
    lgd2.addLine(("cold", "blue", "-", 1))
    lgd2.addLine(("hot", "red", "-", 1))
    # lgd2.addLine(("ISM", "brown", "-", 1))        
    lgd2.addLine(("PhEW", "green", "-", 1))            
    lgd2.draw()
    # lgd3 = legend.legend(axs[5])
    # lgd3.loc = "upper right"
    # lgd3.addLine(("all wind", "orange", "-", 1))        
    # lgd3.addLine(("PhEW", "green", "-", 1))
    # lgd3.addLine(("cold wind", "cyan", "-", 1))
    # lgd3.addLine(("hot wind", "magenta", "-", 1))
    # lgd3.draw()

    lgd3 = legend.legend(axs[1])
    lgd3.loc = "upper right"
    lgd3.addLine(("Averaged", "black", "-", 1))        
    lgd3.addLine(("cold metals", "blue", "-", 1))
    lgd3.addLine(("hot metals", "red", "-", 1))
    # lgd3.addLine(("wind metals", "green", "-", 1))
    lgd3.draw()

plt.savefig(DIRS['FIGURE']+"tmp.pdf")
plt.show()

print ("Done.")
