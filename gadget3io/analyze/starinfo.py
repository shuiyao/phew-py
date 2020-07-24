# Compile Information for star particles.
# Find Mgal and Mvir and specify if center

from cosmology import acosmic, tcosmic
import ioformat
from numpy import genfromtxt, linspace, array, sort, insert
from numpy import where, inf, log10, isinf
import matplotlib.pyplot as plt
import config_mpl

print "compiled."
DRAW_FRAME = False

REDSHIFT = 1.0
if(REDSHIFT == 1.0):
    snapstr = "078"
if(REDSHIFT == 0.0):
    snapstr = "108"

TMAX_CUT = 5.5

def find_fnames(modelname):
    # snapstr = "058"
    galname = "/proj/shuiyao/"+modelname+"/"+"gal_z"+snapstr+".stat"
    soname = "/proj/shuiyao/"+modelname+"/"+"so_z"+snapstr+".sovcirc"
    finname = "/scratch/shuiyao/scidata/gadget3io/"+modelname+"/"+modelname+"_"+snapstr+".stars"
    foutname = "/scratch/shuiyao/scidata/gadget3io/"+modelname+"/"+modelname+"_"+snapstr+".starhost"
    return galname, soname, finname, foutname

HUBBLEPARAM = 0.7
OMEGAB = 0.045
OMEGAM = 0.30

def write_starhost_file(modelname, unit_m):
    # modelname = "l50n288-phew-m5"
    # unit_m = 3.4e6 * 8.
    galname, soname, finname, foutname = find_fnames(modelname)
    stars = genfromtxt(finname, names=True)
    mstar = ioformat.rcol(galname, [4], linestart=0)
    msub = ioformat.rcol(soname, [6], linestart=1)
    mstar = log10(array(mstar) * unit_m * 1.e10 / 0.7)
    msub = log10(array(msub) / 0.7)
    stars['Mass'] *= unit_m * 1.e10 / 0.7
    print "Start Writing"
    fout = open(foutname, "w")
    fout.write("#Mass Mstar Msub Tmax Flag\n")
    for s in stars:
        gidx, hidx = (int)(s['GID']) - 1, (int)(s['HID']) - 1
        if(gidx < 0 or hidx < 0): continue
        ms, mh = mstar[gidx], msub[hidx]
        if(gidx == hidx): flag = 1
        else: flag = 0
        line = "%7.5e %6.3f %6.3f %5.3f %d\n" % \
               (s['Mass'], ms, mh, s['Tmax'], flag)
        fout.write(line)
    fout.close()

models = ['l25n288-phew-m5', 'l25n288-phew-m5']
# models = ['l25n288-phewoff-fw', 'l50n288-phewoff']
clrs = ["red", "orange"]
lstyles = ["-", "--"]
step = 100

def info_starhosts():
    for mi, modelname in enumerate(models):
        foutname = "/scratch/shuiyao/scidata/gadget3io/"+modelname+"/"+modelname+"_"+snapstr+".starhost"
        print foutname
        stars = genfromtxt(foutname, names=True)
        print "Number of Stars = %8d, total Mass = %5.2f" % (len(stars), log10(sum(stars['Mass'])))
        stars = stars[stars['Flag'] == 1]
        print "Number of Central Stars = %8d, total Mass = %5.2f" % (len(stars), log10(sum(stars['Mass'])))
        print "Halo Mass Range: ", min(stars['Msub']), max(stars['Msub'])

def draw_mvbins():
    fig, ax = plt.subplots(1,1, figsize=(8,7))
    for mi, modelname in enumerate(models):
        foutname = "/scratch/shuiyao/scidata/gadget3io/"+modelname+"/"+modelname+"_"+snapstr+".starhost"
        print "Building MBins for ", foutname
        mvbins = build_msub_bins(foutname)
        mvbins.draw(ax, lstyle=lstyles[mi])
    from pltastro import legend
    lgd1 = legend.legend(ax)
    lgd1.loc = "upper left"
    lgd1.addLine((models[0], "black", lstyles[0], 1))
    lgd1.addLine((models[1], "black", lstyles[1], 1))
    lgd1.draw()

        # tab = genfromtxt(foutname, names=True)
        # tab = tab[tab['Flag'] == 1.]
        # plt.plot(tab['Msub'][::step], tab['Tmax'][::step], ".", color=clrs[mi], markersize=6, alpha=0.2)
        # plt.plot(tab['Msub'][::step], tab['Mstar'][::step], ".", color=clrs[mi], markersize=6, alpha=0.2)
    ax.set_title("z = "+str(REDSHIFT)[:3])
    ax.set_xlabel(r"$\log(M_\mathrm{vir}/M_\odot)$")
    ax.set_ylabel("Mass Fraction")
    ax.set_xlim(10.5, 14.0)
    ax.set_ylim(0.0, 1.1)
    plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
    
    plt.show()
    plt.close()

class binned_by_virial_mass():
    def __init__(self, nbins=40, mmin=10.0, mmax=14.0):
        self.nbins = nbins
        self.binmin = mmin
        self.binlen = (mmax - mmin) / (float)(nbins)
        self.nodes = linspace(mmin, mmax, nbins+1)
        self.mid = 0.5 * (self.nodes[:-1] + self.nodes[1:])
        self.cold = array([0.0]*self.nbins)
        self.hot = array([0.0]*self.nbins)
        self.total = array([0.0]*self.nbins)
    def find_idx_for_bin(self, x):
        bidx = (int)((x - self.binmin) / self.binlen)
        if(bidx < 0): bidx = 0
        if(bidx > self.nbins-1): bidx = self.nbins - 1
        return bidx
    def draw(self, ax, lstyle="-"):
        ax.plot(self.mid, self.cold/self.total, linestyle=lstyle, color="blue")
        ax.plot(self.mid, self.hot/self.total, linestyle=lstyle, color="red")
        ax.plot(self.mid, self.total/self.total, linestyle=lstyle, color="black")

def build_msub_bins(foutname):
    mvbins = binned_by_virial_mass(nbins=40, mmin=10.0, mmax=14.0)
    stars = genfromtxt(foutname, names=True)
    stars = stars[stars['Flag'] == 1.]
    for si, s in enumerate(stars): # loop over all stars
        msub = s['Msub']
        bidx = mvbins.find_idx_for_bin(msub)
        mvbins.total[bidx] += s['Mass']    
        if(s['Tmax'] > TMAX_CUT):
            mvbins.hot[bidx] += s['Mass']
        else:
            mvbins.cold[bidx] += s['Mass']
    return mvbins

