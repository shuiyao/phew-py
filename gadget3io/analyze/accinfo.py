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
import matplotlib.pyplot as plt
import config_mpl

Mgas_orig = 1.167e7 # Hires
Mgas_orig = 9.337e7 # Hires

print "compiled."
DRAW_FRAME = True
TMAX_CUT = 5.5
SHOW_ALLWINDS = True
RENEW = True

# modelnames = ["p50n288sw", "l50n288-phewoff", "p50n288beta2"] # mi
# lgds = ["RefSlow", "Ref", "Ref$\sigma$3"]
lstyles = ["--", "-", ":"]
# modelnames = ['l25n288-phew-m4', 'l25n288-phew-m5']
# lgds = ["PhEW-L25-Mc4", "PhEW-L25-Mc5"]
#modelnames = ['l50n288-phew-m4', 'l50n288-phew-m5']
#lgds = ["PhEW-L50-Mc4", "PhEW-L50-Mc5"]
# modelnames = ['p50n288fiducial', 'l50n288-phewoff']
# lgds = ["SPH", "GIZMO"]
# modelnames = ['l50n288-phewoff', 'l25n288-phewoff-fw']
# lgds = ["GIZMO", "GIZMO-Hres"]
# modelnames = ['l50n288-phew-m5', 'l25n288-phew-m5']
# lgds = ["PhEW-L50-Mc5", "PhEW-L25-Mc5"]
# modelnames = ['l25n144-phew-m5', 'l25n288-phew-m5']
# lgds = ["PhEW,25/144", "PhEW,25/288"]
modelnames = ['l25n288-phew-m5', 'l25n288-phew-m5-spl']
lgds = ["PhEW,25/288", "PhEW,25/288,Split"]
# modelnames = ['l25n144-phew-rcloud', 'l50n288-phew-m5']
# lgds = ["PhEW,25/144", "PhEW,50/288"]
# modelnames = ['l25n144-phewoff', 'l25n288-phewoff-fw']
# lgds = ["PhEW-L25N144-Mc5", "PhEW-L50-Mc5"]
REDSHIFT = 0.25
if(REDSHIFT == 0.0):
    zstr = "108"
if(REDSHIFT == 2.0):
    zstr = "058"
if(REDSHIFT == 1.0):
    zstr = "078"
if(REDSHIFT == 0.25):
    zstr = "098"    

def find_fnames(modelname):
    snapstr = zstr
    sfrinfoname = "/scratch/shuiyao/scidata/gadget3io/"+modelname+"/"+modelname+"_"+snapstr+".starinfo"
    soname = "/proj/shuiyao/"+modelname+"/"+"so_z"+snapstr+".sovcirc"
    if(TMAX_CUT == 5.5):
        foutname = "/scratch/shuiyao/scidata/gadget3io/"+modelname+"/"+modelname+"_"+snapstr+".accinfo"
    if(TMAX_CUT == 6.0):        
        foutname = "/scratch/shuiyao/scidata/gadget3io/"+modelname+"/"+modelname+"_"+snapstr+"_t6.accinfo"
    return sfrinfoname, soname, foutname

FIGNAME = "sfhistoryz2.pdf"
color_coldw = "cyan"
color_hotw = "magenta"
color_allw = 'orange'
# color_mixw = "cyan"
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
        self.wmix = array([0.0]*self.nbins)        
        self.other = array([0.0]*self.nbins)
        self.total = array([0.0]*self.nbins)
        self.totalmass = array([0.0]*self.nbins)        
    def find_idx_for_bin(self, x):
        bidx = (int)((x - self.binmin) / self.binlen)
        if(bidx < 0): bidx = 0
        if(bidx > self.nbins-1): bidx = self.nbins - 1
        return bidx
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
    def write(self, fname):
        fout = open(fname, "w")
        fout.write("#Mvir Mcold Mhot Mwcold Mwhot Mtot Mtotalmass\n")
        for i in range(len(self.mid)):
            line = "%g %g %g %g %g %g %g\n" % \
                   (self.mid[i], self.cold[i], self.hot[i], self.wcold[i], \
                    self.whot[i], self.total[i], self.totalmass[i])
            fout.write(line)
        fout.close()
    def draw(self, ax, lstyle="-", fraction=False):
        if(fraction == False):
            ax.plot(self.mid, log10(self.cold/self.totalmass), linestyle=lstyle, color="blue")
            ax.plot(self.mid, log10(self.hot/self.totalmass), linestyle=lstyle, color="red")
            ax.plot(self.mid, log10(self.wcold/self.totalmass), linestyle=lstyle, color=color_coldw)
            ax.plot(self.mid, log10(self.whot/self.totalmass), linestyle=lstyle, color=color_hotw)
            if(SHOW_ALLWINDS):
                ax.plot(self.mid, log10((self.whot+self.wcold)/self.totalmass), linestyle=lstyle, color=color_allw)                
            ax.plot(self.mid, log10((self.total-self.other)/self.totalmass), linestyle=lstyle, color="lightgrey")
            ax.plot(self.mid, log10(self.total/self.totalmass), linestyle=lstyle, color="black")
        else:
            ax.plot(self.mid, self.cold/self.total, linestyle=lstyle, color="blue")
            ax.plot(self.mid, self.hot/self.total, linestyle=lstyle, color="red")
            ax.plot(self.mid, self.wcold/self.total, linestyle=lstyle, color=color_coldw)
            ax.plot(self.mid, self.whot/self.total, linestyle=lstyle, color=color_hotw)
            if(SHOW_ALLWINDS):
                ax.plot(self.mid, (self.whot+self.wcold)/self.total, linestyle=lstyle, color=color_allw)                
            ax.plot(self.mid, 1.0-self.other/self.total, linestyle=lstyle, color="lightgrey")
            ax.plot(self.mid, self.total/self.total, linestyle=lstyle, color="black")

def read_starinfo(fname, fformat=1):
    if(fformat == 1):
        stars = genfromtxt(fname, dtype='f8,f8,f8,f8,f8,f8,f8,f8,f8,i8,i8', names=True)
    else:
        stars = genfromtxt(fname, dtype='f8,f8,f8,f8,f8,f8,f8,i8,i8', names=True)
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
            if(s['Tmax'] > TMAX_CUT):
                mvbins.hot[bidx] += s['Mass'] - s['WindMass']
                mvbins.whot[bidx] += s['WindMass']
            else:
                mvbins.cold[bidx] += s['Mass'] - s['WindMass']
                mvbins.wcold[bidx] += s['WindMass']
        elif(s['a_last'] < 0): # wind
            if(abs(s['Tmax']) > TMAX_CUT): mvbins.whot[bidx] += s['Mass']
            else: mvbins.wcold[bidx] += s['Mass']
        else: mvbins.other[bidx] += s['Mass']
        if(flag[hidx] == 0):
            flag[hidx] = 1
            mvbins.totalmass[bidx] += mvir
    mvbins.totalmass *= (OMEGAB / OMEGAM) * 2.0e33
    return mvbins

def load_central_stars(fname, fformat):
    stars = read_starinfo(fname, fformat)
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
    mtot = mtot * 2.e33    
    return mtot

if(DRAW_FRAME):
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

def draw_frame():
    from matplotlib.lines import Line2D
    from pltastro import legend
    lgd1 = legend.legend(axs[0])
    lgd1.loc = "lower left"
    # lgd1.addLine((lgds[2], "black", ":", 1))
    lgd1.addLine((lgds[0], "black", "--", 1))
    lgd1.addLine((lgds[1], "black", "-", 1))
    lgd1.draw()
    lgd2 = legend.legend(axs[0])
    lgd2.loc = "lower right"
    lgd2.addLine(("cold", "blue", "-", 1))
    lgd2.addLine(("hot", "red", "-", 1))
    lgd2.addLine(("cold wind", "cyan", "-", 1))
    lgd2.addLine(("hot wind", "magenta", "-", 1))
    if(SHOW_ALLWINDS):
        lgd2.addLine(("all wind", color_allw, "-", 1))    
    lgd2.draw()

    FBEHROOZI = "/scratch/shuiyao/sci/REFERENCES/behroozi13/smmr/"
    if(REDSHIFT == 2.0):
        snapnum, zstr_behroozi = "058", "2.00"
    if(REDSHIFT == 0.0):
        snapnum, zstr_behroozi = "108", "0.10"
    if(REDSHIFT == 1.0):
        snapnum, zstr_behroozi = "078", "1.00"
    if(REDSHIFT == 0.25):
        snapnum, zstr_behroozi = "098", "0.10"
    mh, ms, err1, err2 = ioformat.rcol(FBEHROOZI+"c_smmr_z"+zstr_behroozi+"_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat", [0,1,2,3], linestart=1)
    ms = array(ms) - log10(0.156)
    axs[0].plot(mh, ms, "o", color="green", markersize=6)
    axs[0].errorbar(mh, ms, yerr=[err1, err2], color="green")

    axs[0].set_title("z = "+str(REDSHIFT)[:3], fontsize=16)

    axs[0].plot([10.9, 10.9], [-4.0, -0.5], "k--")
    axs[0].plot([10.0, 10.0], [-4.0, -0.5], "k-")
    axs[1].plot([10.9, 10.9], [0.0, 1.1], "k--")
    axs[1].plot([10.0, 10.0], [0.0, 1.1], "k-")

    plt.savefig(FIGNAME)
    plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
    plt.show()
    
fformats = [0, 1]
for mi in range(2):
    sfrinfoname, soname, foutname = find_fnames(modelnames[mi])
    if(not os.path.exists(foutname) and RENEW == False):
        stars = load_central_stars(sfrinfoname, fformat=fformats[mi])
        mvbins = build_mvir_bins(stars, soname) # totalmass estimated here
        mvbins.write(foutname)
    # OBSOLETE: mvbins.totalmass = find_total_mass(stars, soname) 
    else:
        mvbins = binned_by_virial_mass(nbins=40, mmin=10.0, mmax=14.0)
        mvbins.read(foutname)
        mvbins.draw(axs[0], lstyle=lstyles[mi], fraction=False)
        mvbins.draw(axs[1], lstyle=lstyles[mi], fraction=True)
if(DRAW_FRAME == True):
    draw_frame()

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
