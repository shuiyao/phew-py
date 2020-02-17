# Compile Information for star particles.
# Including their:
# Mstar, Mvir, aform, aacc, Tmax, alast

# Sfrinfo:
# a, ID, LastSFTime, Tmax, ...
# LastSFTime =
#   - 0: First Time SF, pristine gas accretion
#   - +: non-SF -> SF
#   - -: Wind -> SF
# Tmax is set to 0 at launch.
from cosmology import acosmic, tcosmic
import ioformat
from numpy import genfromtxt, linspace, array, sort, insert
from numpy import where, inf, log10, isinf
import config_mpl
import matplotlib.pyplot as plt
import config_mpl

print "compiled."
def find_fnames(modelname):
    snapstr = "058"
    sfrinfoname = "/scratch/shuiyao/scidata/gadget3io/"+modelname+"/"+modelname+"_"+snapstr+".starinfo"
    soname = "/scratch/shuiyao/data/"+modelname+"/"+"so_z"+snapstr+".sovcirc"    
    return sfrinfoname, soname

FIGNAME = "sfhistoryz2.pdf"
color_coldw = "cyan"
color_hotw = "magenta"
HUBBLEPARAM = 0.7
OMEGAB = 0.045
OMEGAM = 0.30

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
        self.other = array([0.0]*self.nbins)
        self.total = array([0.0]*self.nbins)
        self.totalmass = array([0.0]*self.nbins)        
    def find_idx_for_bin(self, x):
        bidx = (int)((x - self.binmin) / self.binlen)
        if(bidx < 0): bidx = 0
        if(bidx > self.nbins-1): bidx = self.nbins - 1
        return bidx
    def draw(self, ax, lstyle="-", fraction=False):
        if(fraction == False):
            ax.plot(self.mid, log10(self.cold/self.totalmass), linestyle=lstyle, color="blue")
            ax.plot(self.mid, log10(self.hot/self.totalmass), linestyle=lstyle, color="red")
            ax.plot(self.mid, log10(self.wcold/self.totalmass), linestyle=lstyle, color=color_coldw)
            ax.plot(self.mid, log10(self.whot/self.totalmass), linestyle=lstyle, color=color_hotw)
            ax.plot(self.mid, log10((self.total-self.other)/self.totalmass), linestyle=lstyle, color="lightgrey")
            ax.plot(self.mid, log10(self.total/self.totalmass), linestyle=lstyle, color="black")
        else:
            ax.plot(self.mid, self.cold/self.total, linestyle=lstyle, color="blue")
            ax.plot(self.mid, self.hot/self.total, linestyle=lstyle, color="red")
            ax.plot(self.mid, self.wcold/self.total, linestyle=lstyle, color=color_coldw)
            ax.plot(self.mid, self.whot/self.total, linestyle=lstyle, color=color_hotw)
            ax.plot(self.mid, 1.0-self.other/self.total, linestyle=lstyle, color="lightgrey")
            ax.plot(self.mid, self.total/self.total, linestyle=lstyle, color="black")

# class binned_by_aform():
#     def __init__(self, nbins=40, amin=0.0, amax=1.0):
#         self.nbins = nbins
#         self.binmin = amin
#         self.binlen = (amax - amin) / (float)(nbins)
#         self.nodes = linspace(amin, amax, nbins+1)
#         self.mid = 0.5 * (self.nodes[:-1] + self.nodes[1:])
#         self.cold = array([0.0]*self.nbins)
#         self.hot = array([0.0]*self.nbins)
#         self.wcold = array([0.0]*self.nbins)
#         self.whot = array([0.0]*self.nbins)
#         self.other = array([0.0]*self.nbins)
#         self.total = array([0.0]*self.nbins)
#         self.totalmass = 0.0
#     def find_idx_for_bin(self, x):
#         bidx = (int)((x - self.binmin) / self.binlen)
#         if(bidx < 0): bidx = 0
#         if(bidx > self.nbins-1): bidx = self.nbins - 1
#         return bidx
#     def draw(self, ax, lstyle="-", fraction=False):
#         if(fraction == False):
#             ax.plot(self.mid, log10(self.cold/self.totalmass), linestyle=lstyle, color="blue")
#             ax.plot(self.mid, log10(self.hot/self.totalmass), linestyle=lstyle, color="red")
#             ax.plot(self.mid, log10(self.wcold/self.totalmass), linestyle=lstyle, color=color_coldw)
#             ax.plot(self.mid, log10(self.whot/self.totalmass), linestyle=lstyle, color=color_hotw)
#             ax.plot(self.mid, log10((self.total-self.other)/self.totalmass), linestyle=lstyle, color="lightgrey")
#             ax.plot(self.mid, log10(self.total/self.totalmass), linestyle=lstyle, color="black")
#         else:
#             ax.plot(self.mid, self.cold/self.total, linestyle=lstyle, color="blue")
#             ax.plot(self.mid, self.hot/self.total, linestyle=lstyle, color="red")
#             ax.plot(self.mid, self.wcold/self.total, linestyle=lstyle, color=color_coldw)
#             ax.plot(self.mid, self.whot/self.total, linestyle=lstyle, color=color_hotw)
#             ax.plot(self.mid, 1.0-self.other/self.total, linestyle=lstyle, color="lightgrey")
#             ax.plot(self.mid, self.total/self.total, linestyle=lstyle, color="black")

def read_starinfo(fname):
    stars = genfromtxt(fname, dtype='f8,f8,f8,f8,f8,f8,i8,i8', names=True)    
    # a_form a_acc a_last Mass Tmax Mstar Mvir
    # if Mvir < 0, it's a satellite galaxy
    return stars

def build_mvir_bins(stars, soname):
    msub = ioformat.rcol(soname, [6], linestart=1)
    flag = [0] * len(msub)
    mvbins = binned_by_virial_mass(nbins=40, mmin=10.0, mmax=14.0)
    for s in stars: # loop over all stars
        hidx = s['HID']-1
        mvir = msub[hidx] / HUBBLEPARAM
        if(isinf(log10(mvir))): continue
        bidx = mvbins.find_idx_for_bin(log10(mvir))
        mvbins.total[bidx] += s['Mass']    
        if(s['a_last'] == 0): # primordial
            if(s['Tmax'] > 5.5): mvbins.hot[bidx] += s['Mass']
            else: mvbins.cold[bidx] += s['Mass']        
        elif(s['a_last'] < 0): # wind
            if(s['Tmax'] < -5.5): mvbins.whot[bidx] += s['Mass']
            else: mvbins.wcold[bidx] += s['Mass']
        else: mvbins.other[bidx] += s['Mass']
        if(flag[hidx] == 0):
            flag[hidx] = 1
            mvbins.totalmass[bidx] += mvir
    mvbins.totalmass *= (OMEGAB / OMEGAM)
    return mvbins

def load_central_stars(fname):
    stars = read_starinfo(fname)
    nstars = len(stars)
    stars = stars[stars['GID'] == stars['HID']] # Only central galaxies    
    print "Stars from central galaxies: %d/%d (%4.1f%%)" % \
        (len(stars), nstars, (float)(len(stars)) / (float)(nstars) * 100.)
    return stars

def find_total_mass(stars, soname):
    mtot = 0.0
    msub = ioformat.rcol(soname, [6], linestart=1)
    ncount = 0
    for s in stars:
        m = msub[s['HID'] - 1]
        if(m > 0):
            mtot += m
            ncount += 1
            msub[s['HID'] - 1] = 0.0 # Once added, reset to 0 to prevent repeated counting
    mtot = mtot / HUBBLEPARAM
    print "Total amount of haloes found: %d." % (ncount)
    print "  - Average mass (log): %5.2f." % log10((mtot / float(ncount)))
    mtot = mtot * (OMEGAB / OMEGAM)
    return mtot

from pltastro import *
import pltastro
frm = pltastro.frame.multi(2, 1)
frm.panels.ylabels[0] = r"$\log(\frac{M_*}{(\Omega_b/\Omega_m)M_{h}})$"
frm.panels.ylabels[1] = r"$f_{acc}$"
frm.panels.set_xlabels(r"$\log\ M_{vir}$")
frm.panels.set_xlims(10., 14.)
frm.panels.ylims[1] = (0.0, 1.1)
frm.panels.yticks[1] = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
frm.panels.ylims[0] = (-4.0, -0.5)
frm.panels.yticks[0] = [-4.0, -3.0, -2.0, -1.0]
frm.panels.ytickformat = ['%3.1f', "%3.1f"]
frm.panels.set_xlims(11.0, 13.5)
frm.params.figsize = (6, 7)
frm.params.height_ratios = [2, 1]
frm.params.left = 0.20
fig, axs = draw(frm)

modelnames = ["p50n288sw", "p50n288fiducial", "p50n288beta2"] # mi
lgds = ["RefSlow", "Ref", "Ref$\sigma$3"]
lstyles = ["--", "-", ":"]

for mi in range(3):
    sfrinfoname, soname = find_fnames(modelnames[mi])
    stars = load_central_stars(sfrinfoname)
    mvbins = build_mvir_bins(stars, soname) # totalmass estimated here
    # mvbins.totalmass = find_total_mass(stars, soname)    
    mvbins.draw(axs[0], lstyle=lstyles[mi], fraction=False)
    mvbins.draw(axs[1], lstyle=lstyles[mi], fraction=True)

from matplotlib.lines import Line2D
from pltastro import legend
lgd1 = legend.legend(axs[0])
lgd1.loc = "lower left"
lgd1.addLine((lgds[2], "black", ":", 1))
lgd1.addLine((lgds[0], "black", "--", 1))
lgd1.addLine((lgds[1], "black", "-", 1))
lgd1.draw()
lgd2 = legend.legend(axs[0])
lgd2.loc = "lower right"
lgd2.addLine(("cold", "blue", "-", 1))
lgd2.addLine(("hot", "red", "-", 1))
lgd2.addLine(("cold wind", "cyan", "-", 1))
lgd2.addLine(("hot wind", "magenta", "-", 1))
lgd2.draw()
# lgds = [Line2D([0], [0], color="black", linestyle=":", label=lgds[2]),\
#         Line2D([0], [0], color="black", linestyle="--", label=lgds[0]),\
#         Line2D([0], [0], color="black", linestyle="-", label=lgds[1])]
# l1 = axs[0].legend(handles=lgds, loc="lower left")
# axs[0].add_artist(l1)
# lgds = [Line2D([0], [0], color="blue", linestyle="-", label="cold"),\
#         Line2D([0], [0], color="red", linestyle="-", label="hot"),\
#         Line2D([0], [0], color="cyan", linestyle="-", label="cold wind"),\
#         Line2D([0], [0], color="magenta", linestyle="-", label="hot wind")]
# axs[0].legend(handles=lgds, loc="lower right")

# FMOSTER = "/scratch/shuiyao/sci/REFERENCES/moster_2013_2018/"
# zstr_moster = "2"
# mh, ms = ioformat.rcol(FMOSTER+"moster18_z"+zstr_moster+".dat", [0,1])
# ms = array(ms) - (array(mh) + log10(0.156))
# axs[0].plot(mh, ms, "-", color="purple", linewidth=2)

FBEHROOZI = "/scratch/shuiyao/sci/REFERENCES/behroozi13/smmr/"
snapnum, zstr_behroozi = "058", "2.00"
mh, ms, err1, err2 = ioformat.rcol(FBEHROOZI+"c_smmr_z"+zstr_behroozi+"_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat", [0,1,2,3], linestart=1)
ms = array(ms) - log10(0.156)
axs[0].plot(mh, ms, "o", color="green", markersize=6)
axs[0].errorbar(mh, ms, yerr=[err1, err2], color="green")

axs[0].set_title("z = 2", fontsize=16)
plt.savefig(FIGNAME)

def show_wind_temperatures():
    winds = stars[stars['a_last'] < 0]
    # plt.plot(winds['a_form'] - abs(winds['a_last']), winds['Tmax'], "b.", alpha=0.05)
    plt.plot(stars['Mvir'], abs(stars['Tmax']), "b.", alpha=0.01, markersize=2)
    plt.plot(winds['Mvir'], abs(winds['Tmax']), "r.", alpha=0.01, markersize=2)
    import cosmology
    from astroconst import pc, ac
    mhalos = linspace(10., 14., 41)
    Tvirs = log10(cosmology.Tvir(10.**mhalos * ac.msolar, z = 2.0))
    plt.plot(mhalos, Tvirs, "k--")
    plt.show()
