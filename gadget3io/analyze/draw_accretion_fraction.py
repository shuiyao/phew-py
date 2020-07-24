# Making a specific plot for Neal
from cosmology import acosmic, tcosmic
from numpy import genfromtxt, log10, linspace, array
import matplotlib.pyplot as plt
import config_mpl

from pltastro import *
import pltastro

color_coldw = "cyan"
color_hotw = "magenta"
HUBBLEPARAM = 0.7
OMEGAB = 0.045
OMEGAM = 0.30

fig, ax = plt.subplots(1, 1, figsize=(7,6))
axs = [ax]

def find_fnames(modelname):
    snapstr = "108"
    sfrinfoname = "/scratch/shuiyao/scidata/gadget3io/"+modelname+"/"+modelname+"_"+snapstr+".starinfo"
    soname = "/proj/shuiyao/"+modelname+"/"+"so_z"+snapstr+".sovcirc"
    foutname = "/scratch/shuiyao/scidata/gadget3io/"+modelname+"/"+modelname+"_"+snapstr+".accinfo"
    return sfrinfoname, soname, foutname

class binned_by_virial_mass():
    def __init__(self, nbins=40, mmin=10.0, mmax=14.0):
        self.nbins = nbins
        self.binmin = mmin
        self.binlen = (mmax - mmin) / (float)(nbins)
        self.nodes = linspace(mmin, mmax, nbins+1)
        self.mid = 0.5 * (self.nodes[:-1] + self.nodes[1:])
        self.cold = array([0.0]*self.nbins)
        self.hot = array([0.0]*self.nbins)
        self.wcold = array([0.0]*self.nbins)
        self.whot = array([0.0]*self.nbins)
        self.wmix = array([0.0]*self.nbins)        
        self.other = array([0.0]*self.nbins)
        self.total = array([0.0]*self.nbins)
        self.totalmass = array([0.0]*self.nbins)        
    def read(self, fname):
        tab = genfromtxt(fname, names=True)
        for i in range(self.nbins):
            self.mid[i] = tab['Mvir'][i]
            self.cold[i] = tab['Mcold'][i]
            self.hot[i] = tab['Mhot'][i]
            self.wcold[i] = tab['Mwcold'][i]
            self.whot[i] = tab['Mwhot'][i]
            self.total[i] = tab['Mtot'][i]
            self.totalmass[i] = tab['Mtotalmass'][i]
            self.other[i] = self.total[i] - self.cold[i] - self.hot[i] - self.wcold[i] - self.whot[i] - self.wmix[i]
    def draw(self, ax, lstyle="-"):
        ax.plot(self.mid, self.cold/self.total, linestyle=lstyle, color="blue")
        ax.plot(self.mid, self.hot/self.total, linestyle=lstyle, color="red")
        ax.plot(self.mid, self.wcold/self.total, linestyle=lstyle, color=color_coldw)
        ax.plot(self.mid, self.whot/self.total, linestyle=lstyle, color=color_hotw)
        # ax.plot(self.mid, 1.0-self.other/self.total, linestyle=lstyle, color="lightgrey")
        # ax.plot(self.mid, self.total/self.total, linestyle=lstyle, color="black")

modelnames = ['l25n288-phew-m5', 'l50n288-phew-m5']
lgds = ["PhEW-L25-Mc5", "PhEW-L50-Mc5"]
lstyles = ["--", "-", ":"]

for mi in range(2):
    sfrinfoname, soname, foutname = find_fnames(modelnames[mi])
    mvbins = binned_by_virial_mass(nbins=40, mmin=10.0, mmax=14.0)
    mvbins.read(foutname)
    mvbins.draw(axs[0], lstyle=lstyles[mi])

from matplotlib.lines import Line2D
from pltastro import legend
lgd1 = legend.legend(axs[0])
lgd1.loc = "upper left"
# lgd1.addLine((lgds[2], "black", ":", 1))
lgd1.addLine((lgds[0], "black", "--", 1))
lgd1.addLine((lgds[1], "black", "-", 1))
lgd1.draw()
lgd2 = legend.legend(axs[0])
lgd2.loc = "upper right"
lgd2.addLine(("cold", "blue", "-", 1))
lgd2.addLine(("hot", "red", "-", 1))
lgd2.addLine(("cold wind", "cyan", "-", 1))
lgd2.addLine(("hot wind", "magenta", "-", 1))
lgd2.draw()

axs[0].set_title("z = 0", fontsize=16)
axs[0].set_xlim(11.0, 13.5)
axs[0].set_ylim(0.0, 1.1)
axs[0].set_yticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
axs[0].set_ylabel(r"$f_{acc}$")
axs[0].set_xlabel(r"$\log\ M_{vir}$")

plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
plt.show()

