from mymod import *
from numpy import genfromtxt

hparam = 0.7
unit_m = 1.e10 / hparam

nbins = 30
redges = linspace(0., 1., nbins+1)
dr = 1. / nbins

ZSOLAR = log10(0.0122)

#MODE = "Metal"
MODE = "Mass"
SUB_MODE = "Tmax"
#SUB_MODE = "logT"
RENEW = False


def find_radial_bin(r):
    if(r > 1.0): return -1
    bin_idx = r / dr
    return (int)(bin_idx)

# ---- GIZMO Vs. PhEW
# models = ["l50n288-phewoff", "l50n288-phew-m5"]
# lgds = ["l50n288-phewoff", "l50n288-phew-m5"]
# ---- Mc=1.e4 Vs. Mc=1.e5
# models = ["l50n288-phew-m4", "l50n288-phew-m5"]
# lgds = ["l50n288-phew-m4", "l50n288-phew-m5"]
# ---- 25/288 Vs. 50/288 (PhEW)
# models = ["l50n288-phew-m5", "l25n288-phew-m5"]
# lgds = ["l50n288-phew-m5", "l25n288-phew-m5"]
# ---- 25/288 Vs. 50/288 (GIZMO)
# models = ["l50n288-phewoff", "l25n288-phewoff-fw"]
# lgds = ["l50n288-phewoff", "l25n288-phewoff-fw"]
# ---- 25/144 Vs. 25/288 (PhEW)
models = ["l25n144-phew-m5", "l25n288-phew-m5"]
lgds = ["PhEW,25/144", "PhEW,25/288"]
# models = ["l25n144-phew-m5", "l25n144-phew-m5-spl"]
# lgds = ["PhEW,25/144", "PhEW,25/144,Split"]
lstyles = ["--", "-"]
REDSHIFT = 1.0
if(REDSHIFT == 0.0): zstr = "108"
if(REDSHIFT == 1.0): zstr = "078"
if(REDSHIFT == 2.0): zstr = "058"
if(REDSHIFT == 4.0): zstr = "033"

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
if(MODE == "Mass"):
    panels.ylabels[1] = "Mass Fraction"
if(MODE == "Metal"):
    panels.ylabels[1] = r"$\log(Z/Z_\odot)$"

fig, axs = draw(frm)

captions = [
    r"$11.0 < M_\mathrm{vir} < 11.5$",\
    r"$11.85 < M_\mathrm{vir} < 12.15$",\
    r"$12.85 < M_\mathrm{vir} < 13.15$"\    
]

if(MODE == "Mass"):
    for modeli in range(len(models)):
        for mi, mstr in enumerate(["mh11", "mh12", "mh13"]):
            fname = "/scratch/shuiyao/scidata/gadget3io/"+models[modeli]+"/"+models[modeli]+"_"+zstr+".gas." + mstr
            foutname = "/scratch/shuiyao/scidata/gadget3io/"+models[modeli]+"/"+models[modeli]+"_"+zstr+".mprof." + mstr

            print "Doing: ", fname

            if(os.path.exists(foutname) and RENEW == False):
                print "Existing Mass Profile File."
                tab = genfromtxt(foutname, names=True)
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
            else:
                print "Generating New Mass Profile File."
                tab = genfromtxt(fname, names=True)

                mass = array([0.0] * (nbins + 1))
                mwind = array([0.0] * (nbins + 1))
                mmix = array([0.0] * (nbins + 1))
                mwcold = array([0.0] * (nbins + 1))
                mwhot = array([0.0] * (nbins + 1))                
                mcold = array([0.0] * (nbins + 1))
                mhot = array([0.0] * (nbins + 1))
                mism = array([0.0] * (nbins + 1))
                mmaxhot = array([0.0] * (nbins + 1))
                mmaxcold = array([0.0] * (nbins + 1))
                mmaxwhot = array([0.0] * (nbins + 1))
                mmaxwcold = array([0.0] * (nbins + 1))                

                for i in range(len(tab)):
                    part = tab[i]
                    bidx = find_radial_bin(part['dr'] / part['Rvir'])
                    mass[bidx] += part['Mass'] # Total Mass
                    mmix[bidx] += part['WMass'] # Total Mixed Mass
                    if(part['SfFlag'] == 1):
                        mism[bidx] += part['Mass']
                    else:
                        if(part['Mc'] > 0): # A PhEW Particle
                            mwind[bidx] += part['Mass']
                            mmix[bidx] -= part['WMass']
                        if(part['Mc'] == 0): # A normal gas particle
                            if(part['logT'] > 5.5):
                                mhot[bidx] += part['Mass'] - part['WMass']
                                mwhot[bidx] += part['WMass']
                            else:
                                mcold[bidx] += part['Mass'] - part['WMass']
                                mwcold[bidx] += part['WMass']
                            if(part['Tmax'] > 5.5):
                                mmaxhot[bidx] += part['Mass'] - part['WMass']
                                mmaxwhot[bidx] += part['WMass']
                            else:
                                mmaxcold[bidx] += part['Mass'] - part['WMass']
                                mmaxwcold[bidx] += part['WMass']
                        if(part['Mc'] < 0): # Tricky, now very rare
                            mwind[bidx] += part['Mass'] - part['WMass']
                fwind = mwind/mass
                fout = open(foutname, "w")
                fout.write("#r Mass Mcold Mhot Mwind Mwcold Mwhot Mism Mmaxcold Mmaxhot Mmaxwcold Mmaxwhot\n")
                for i in range(len(redges)):
                    line = "%5.3f %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e\n" % \
                           (redges[i], mass[i], mcold[i], mhot[i], mwind[i], mwcold[i], mwhot[i], mism[i], mmaxcold[i], mmaxhot[i], mmaxwcold[i], mmaxwhot[i])
                    fout.write(line)
                fout.close()

            axs[mi].plot(redges, mass/mass, color="black", linestyle=lstyles[modeli])
            if(SUB_MODE == "logT"):
                axs[mi].plot(redges, mcold/mass, color="blue", linestyle=lstyles[modeli])
                axs[mi].plot(redges, mhot/mass, color="red", linestyle=lstyles[modeli])
                axs[mi].plot(redges, mwcold/mass, color="cyan", linestyle=lstyles[modeli])
                axs[mi].plot(redges, mwhot/mass, color="magenta", linestyle=lstyles[modeli])
            if(SUB_MODE == "Tmax"):
                axs[mi].plot(redges, mmaxcold/mass, color="blue", linestyle=lstyles[modeli])
                axs[mi].plot(redges, mmaxhot/mass, color="red", linestyle=lstyles[modeli])
                axs[mi].plot(redges, mmaxwcold/mass, color="cyan", linestyle=lstyles[modeli])
                axs[mi].plot(redges, mmaxwhot/mass, color="magenta", linestyle=lstyles[modeli])
            axs[mi].plot(redges, mwind/mass, color="green", linestyle=lstyles[modeli])
            axs[mi].plot(redges, mism/mass, color="brown", linestyle=lstyles[modeli])            
            axs[mi].text(0.01, 0.8, captions[mi], transform=axs[mi].transAxes, fontsize=12)
    axs[0].set_title("Radial Profiles, z="+str(REDSHIFT)[:3])

    from pltastro import legend
    lgd1 = legend.legend(axs[0])
    lgd1.loc = "upper right"
    lgd1.addLine((lgds[0], "black", lstyles[0], 1))
    lgd1.addLine((lgds[1], "black", lstyles[1], 1))
    lgd1.draw()
    lgd2 = legend.legend(axs[2])
    lgd2.loc = "upper right"
    lgd2.addLine(("cold", "blue", "-", 1))
    lgd2.addLine(("hot", "red", "-", 1))
    lgd2.addLine(("cold wind", "cyan", "-", 1))
    lgd2.addLine(("hot wind", "magenta", "-", 1))
    lgd2.addLine(("PhEW", "green", "-", 1))    
    lgd2.addLine(("all wind", "orange", "-", 1))        
    lgd2.draw()


# ---------------- METAL ----------------
if(MODE == "Metal"):
    for modeli in range(len(models)):
        for mi, mstr in enumerate(["mh11", "mh12", "mh13"]):
            fname = "/scratch/shuiyao/scidata/gadget3io/"+models[modeli]+"/"+models[modeli]+"_"+zstr+".gas." + mstr
            print "Doing: ", fname
            # if(modeli == 1 or mi == 1): break;

            tab = genfromtxt(fname, names=True)

            mass = array([0.0] * (nbins + 1))
            mz = array([0.0] * (nbins + 1))
            m = array([0.0] * (nbins + 1))
            mzcold = array([0.0] * (nbins + 1))
            mzhot = array([0.0] * (nbins + 1))                

            for i in range(len(tab)):
                part = tab[i]
                bidx = find_radial_bin(part['dr'] / part['Rvir'])
                mass[bidx] += part['Mass'] # Total Mass
                mass_z = part['Mass'] * part['Z'] # Total Metal Mass
                mz[bidx] += mass_z
                if(part['logT'] > 5.5):
                    mzhot[bidx] += mass_z
                else:
                    mzcold[bidx] += mass_z

            axs[mi].plot(redges, log10(mz/mass)-ZSOLAR, color="black", linestyle=lstyles[modeli])
            axs[mi].plot(redges, log10(mzcold/mass)-ZSOLAR, color="blue", linestyle=lstyles[modeli])
            axs[mi].plot(redges, log10(mzhot/mass)-ZSOLAR, color="red", linestyle=lstyles[modeli])
            axs[mi].text(0.50, 0.8, captions[mi], transform=axs[mi].transAxes, fontsize=12)
    axs[0].set_title("Radial Profiles, z="+str(REDSHIFT)[:3])
    for i in range(3):
        axs[i].set_ylim(-2.0, 0.5)
        axs[i].set_yticks([-2.0, -1.5, -1.0, -0.5, 0.0])

    from pltastro import legend
    lgd = legend.legend(axs[2])
    lgd.loc = "lower right"
    lgd.addLine((lgds[0], "black", lstyles[0], 1))
    lgd.addLine((lgds[1], "black", lstyles[1], 1))
    lgd.draw()

# ---------------- G3-WIND ----------------

# for mi, mstr in enumerate(["mh11", "mh12", "mh13"]):
#     fname = "/scratch/shuiyao/scidata/gadget3io/"+models[1]+"/"+models[1]+"_078.gas." + mstr    

#     tab = genfromtxt(fname, names=True)

#     mass = array([0.0] * (nbins + 1))
#     mwind = array([0.0] * (nbins + 1))
#     mmix = array([0.0] * (nbins + 1))

#     for i in range(len(tab)):
#         part = tab[i]
#         bidx = find_radial_bin(part['dr'] / part['Rvir'])
#         mass[bidx] += part['Mass']
#         if(part['Mc'] < 0):            
#             mwind[bidx] += part['Mass']
#     fwind = mwind/mass

#     axs[mi].plot(redges, mass/mass, "k--")
#     axs[mi].plot(redges, mwind/mass, "g--")
#     # axs[mi].text(0.2, 0.8, captions[mi], transform=axs[mi].transAxes, fontsize=16)

plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
plt.show()

print "Done."
