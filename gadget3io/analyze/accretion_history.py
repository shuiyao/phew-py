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
from astroconst import pc, ac

print "compiled."
def find_fnames(modelname, mstr):
    snapstr = "098"
    sfrinfoname = "/scratch/shuiyao/scidata/gadget3io/"+modelname+"/"+modelname+"_"+snapstr+".starinfo."+mstr
    soname = "/proj/shuiyao/"+modelname+"/"+"so_z"+snapstr+".sovcirc"    
    return sfrinfoname, soname

SEPARATE_COLD_HOT_WINDS = True
SHOW_WINDS_MIXED = False

FIGNAME = "sfhistory.pdf"
FIGNAME = "/scratch/shuiyao/figures/tmp.pdf"
CENTRAL_ONLY = True
color_coldw = "cyan"
color_hotw = "magenta"
color_wind = "green"
color_wmix = "cyan"
color_allw = "orange"
HUBBLEPARAM = 0.7
OMEGAB = 0.045
OMEGAM = 0.3

def arraysum(arr): # Turn histogram to cumulative distribution
    for i in range(len(arr))[1:]:
        arr[i] = arr[i] + arr[i-1]
    return arr

class binned_by_aform():
    def __init__(self, nbins=40, amin=0.0, amax=1.0):
        self.nbins = nbins
        self.binmin = amin
        self.binlen = (amax - amin) / (float)(nbins)
        self.nodes = linspace(amin, amax, nbins+1)
        self.mid = log10(1./(0.5 * (self.nodes[:-1] + self.nodes[1:]))) # Redshift
        self.dz = - 1./self.nodes[1:] + 1./self.nodes[:-1]
        self.dz[0] = 999. - 1./self.nodes[1]  
        self.cold = array([0.0]*self.nbins)
        self.hot = array([0.0]*self.nbins)
        self.wcold = array([0.0]*self.nbins)
        self.whot = array([0.0]*self.nbins)
        self.wind = array([0.0]*self.nbins) # Total winds (excluding mixed)
        self.wmix = array([0.0]*self.nbins) # Total winds mixed 
        self.cold = array([0.0]*self.nbins)        
        self.other = array([0.0]*self.nbins)
        self.total = array([0.0]*self.nbins)
        self.totalmass = 0.0 # Total mass of all haloes involved (will be corrected by Ob/Om)
    def find_idx_for_bin(self, x):
        bidx = (int)((x - self.binmin) / self.binlen)
        if(bidx < 0): bidx = 0
        if(bidx > self.nbins-1): bidx = self.nbins - 1
        return bidx
    def cumulate(self):
        self.cold = arraysum(self.cold)
        self.hot = arraysum(self.hot)
        self.wcold = arraysum(self.wcold)
        self.whot = arraysum(self.whot)
        self.wmix = arraysum(self.wmix)
        self.wind = arraysum(self.wind)                
        self.total = arraysum(self.total)
        self.other = arraysum(self.other)        
    def draw(self, ax, lstyle="-", method="cumulative"):
        if(method == "normed"):
            ax.plot(self.mid, log10(self.cold/self.totalmass), linestyle=lstyle, color="blue")
            ax.plot(self.mid, log10(self.hot/self.totalmass), linestyle=lstyle, color="red")
            if(SEPARATE_COLD_HOT_WINDS):
                ax.plot(self.mid, log10(self.wcold/self.totalmass), linestyle=lstyle, color=color_coldw)
                ax.plot(self.mid, log10(self.whot/self.totalmass), linestyle=lstyle, color=color_hotw)
                ax.plot(self.mid, log10((self.whot+self.wcold)/self.totalmass), linestyle=lstyle, color=color_allw)                
            ax.plot(self.mid, log10(self.wind/self.totalmass), linestyle=lstyle, color=color_wind)
            if(SHOW_WINDS_MIXED):
                ax.plot(self.mid, log10(self.wmix/self.totalmass), linestyle=lstyle, color=color_wmix)
                ax.plot(self.mid, log10((self.total-self.wmix)/self.totalmass), linestyle=lstyle, color="lightgrey")
            ax.plot(self.mid, log10((self.total-self.other)/self.totalmass), linestyle=lstyle, color="grey")
            ax.plot(self.mid, log10(self.total/self.totalmass), linestyle=lstyle, color="black")
        elif(method == "differential"):
            ax.plot(self.mid, log10(self.cold/self.totalmass/self.dz), linestyle=lstyle, color="blue")
            ax.plot(self.mid, log10(self.hot/self.totalmass/self.dz), linestyle=lstyle, color="red")
            if(SEPARATE_COLD_HOT_WINDS):            
                ax.plot(self.mid, log10(self.wcold/self.totalmass/self.dz), linestyle=lstyle, color=color_coldw)
                ax.plot(self.mid, log10(self.whot/self.totalmass/self.dz), linestyle=lstyle, color=color_hotw)
                ax.plot(self.mid, log10((self.whot+self.wcold)/self.totalmass/self.dz), linestyle=lstyle, color=color_allw)                
            ax.plot(self.mid, log10(self.wind/self.totalmass/self.dz), linestyle=lstyle, color=color_wind)
            if(SHOW_WINDS_MIXED):
                ax.plot(self.mid, log10(self.wmix/self.totalmass/self.dz), linestyle=lstyle, color=color_wmix)
            ax.plot(self.mid, log10(self.total/self.totalmass/self.dz), linestyle=lstyle, color="black")
            ax.plot(self.mid, log10(1.0-self.other/self.totalmass/self.dz), linestyle=lstyle, color="grey")
        elif(method == "cumulative"):
            ax.plot(self.mid, self.cold/self.total, linestyle=lstyle, color="blue")
            ax.plot(self.mid, self.hot/self.total, linestyle=lstyle, color="red")
            if(SEPARATE_COLD_HOT_WINDS):            
                ax.plot(self.mid, self.wcold/self.total, linestyle=lstyle, color=color_coldw)
                ax.plot(self.mid, self.whot/self.total, linestyle=lstyle, color=color_hotw)
                ax.plot(self.mid, (self.whot+self.wcold)/self.total, linestyle=lstyle, color=color_allw)                
            ax.plot(self.mid, self.wind/self.total, linestyle=lstyle, color=color_wind)
            if(SHOW_WINDS_MIXED):
                ax.plot(self.mid, self.wmix/self.total, linestyle=lstyle, color=color_wmix)
            ax.plot(self.mid, 1.0-self.other/self.total, linestyle=lstyle, color="grey")
            ax.plot(self.mid, self.total/self.total, linestyle=lstyle, color="black")
            ax.plot(self.mid, self.total/self.total[-1], linestyle=lstyle, color="grey") # The stellar mass growth of all stars

def read_starinfo(fname, fformat="new"):
    if(fformat == "new"):
        stars = genfromtxt(fname, dtype='f8,f8,f8,f8,f8,f8,f8,f8,i8,i8', names=True)
    if(fformat == "old"):
        stars = genfromtxt(fname, dtype='f8,f8,f8,f8,f8,f8,i8,i8', names=True)        
    # New Format: a_form a_acc a_last Mass WindMass Mstar Tmax Z GID HID
    # Old Format: a_form a_acc a_last Mass Tmax GID HID    
    # if Mvir < 0, it's a satellite galaxy
    return stars

def build_abins(stars, cumulative=False):
    abins = binned_by_aform(nbins=40, amin=0.0, amax=1.0)
    for s in stars: # loop over all stars
        bidx = abins.find_idx_for_bin(s['a_form'])
        abins.total[bidx] += s['Mass']
        if(s['a_last'] == 0): # primordial
            if(s['Tmax'] > 5.5): abins.hot[bidx] += s['Mass']
            else: abins.cold[bidx] += s['Mass']        
        elif(s['a_last'] < 0): # wind
            # Caution ~ Now Tmax are all > 0
            abins.wind[bidx] += s['Mass']            
            if(s['Tmax'] > 5.5): abins.whot[bidx] += s['Mass']
            else: abins.wcold[bidx] += s['Mass']
        else: abins.other[bidx] += s['Mass']
    if(cumulative == True):
        abins.cumulate()
    abins.hot *= ac.msolar
    abins.cold *= ac.msolar
    abins.wind *= ac.msolar
    abins.wcold *= ac.msolar
    abins.whot *= ac.msolar
    abins.other *= ac.msolar
    abins.total *= ac.msolar    
    return abins

def build_abins_phew(stars, cumulative=False):
    abins = binned_by_aform(nbins=40, amin=0.0, amax=1.0)
    for s in stars: # loop over all stars
        bidx = abins.find_idx_for_bin(s['a_form'])
        abins.total[bidx] += s['Mass']
        # if(s['Mass'] < s['WindMass']): s['WindMass'] = s['Mass']
        if(s['a_last'] == 0): # primordial
            abins.wmix[bidx] += s['WindMass']
            if(s['Tmax'] > 5.5):
                abins.hot[bidx] += s['Mass'] - s['WindMass']
                # if(s['WindMass'] <= 1.0 * s['Mass']):                                
                abins.whot[bidx] += s['WindMass']
            else:
                abins.cold[bidx] += s['Mass'] - s['WindMass']
                # if(s['WindMass'] <= 1.0 * s['Mass']):                
                abins.wcold[bidx] += s['WindMass']
        elif(s['a_last'] < 0): # wind
            # Caution ~ Now Tmax are all > 0
            # abins.wind[bidx] += s['Mass']
            if(s['Tmax'] > 5.5): abins.whot[bidx] += s['Mass']
            else: abins.wcold[bidx] += s['Mass']
        else:
            abins.other[bidx] += s['Mass']
    if(cumulative == True):
        abins.cumulate()
    return abins

def load_central_stars(fname, fformat="new"):
    stars = read_starinfo(fname, fformat=fformat)
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
    mtot = mtot * (OMEGAB / OMEGAM) * ac.msolar
    return mtot

def figure_sfhistory():
    # Show how stars accumulated to form these galaxies at z = 0
    from pltastro import draw
    import pltastro
    frm = pltastro.frame.multi(3, 3)
    # frm.panels.ylabels[0] = r"$\frac{M(a)}{M(1)}$"
    # frm.panels.ylabels[3] = r"$\frac{M}{M_{tot}}$"
    frm.panels.ylabels[0] = r"$\log(\frac{M_*(<z)}{M_b})$"    
    frm.panels.ylabels[3] = r"$\log(\frac{\Delta M_*(z)}{M_b\Delta_z})$"    
    frm.panels.ylabels[6] = r"$f_{acc}$"
    frm.panels.set_xlabels(r"$z_{sf}$")

    for i in [6, 7, 8]: # bottom panels
        frm.panels.ylims[i] = (0.0, 1.1)
        frm.panels.yticks[i] = [0.2, 0.4, 0.6, 0.8, 1.0]
    for i in [0, 1, 2]: # upper panels
        frm.panels.ylims[i] = (-4.0, 0.0)
        frm.panels.yticks[i] = [-4.0, -3.0, -2.0, -1.0]
    for i in [3, 4, 5]: # middle panels
        frm.panels.ylims[i] = (-5.0, 0.0)
        frm.panels.yticks[i] = [-5.0, -4.0, -3.0, -2.0, -1.0]            

    # frm.panels.set_ylims(0.0, 1.1)
    # frm.panels.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # frm.panels.set_xticks([0.0, 0.33, 0.5, 0.8])
    frm.panels.set_xticks([log10(5.0), log10(3.0), log10(2.0), 0.0])    
    frm.panels.set_xticklabels(["4.0", "2.0", "1.0", "0.0"])    
    frm.panels.set_ytickformat('%3.1f')
    frm.panels.set_xlims(log10(7.0), 0.0)
    frm.params.figsize = (9, 9)
    frm.params.left = 0.12
    fig, axs = draw(frm)

    mstrs = ['mh11', 'mh12', 'mh13'] # fi
    # modelnames = ["p50n288sw", "p50n288fiducial", "p50n288beta2"] # mi

    # modelnames = ["l25n144-gadget3","l25n144-phewoff", "l25n144-phew"] # mi
    # labels = ["G3Winds", "PhEWOFF", "PhEW"]    

    # modelnames = ["l25n144-phew-m5kh100fs10","l25n144-phew-m4kh50fs10", "l25n144-phew-m4kh100fs10"] # mi    
    # labels = ["mc5", "mc4kh50", "mc4kh100"]

    # modelnames = ["l25n144-phewoff","l25n144-phew-m5kh30fs10", "l25n144-phew-mach1"] # mi    
    # labels = ["PhEWOff", "m5kh30fs10", "m5kh30fs10"]

    # modelnames = ["l25n144-phew-m4kh100fs100","l25n144-phew-m4kh50fs10", "l25n144-phew-m4kh100fs10"] # mi    
    # labels = ["mc4fs100", "mc4kh50", "mc4fs10"]

    # modelnames = ["p50n288fiducial","l25n144-gadget3", "l25n144-phewoff"] # mi
    # labels = ["p50n288", "Gadget3", "GIZMO"]

    # modelnames = ["l25n288-phew-m5","l25n144-phew-m5"] # mi
    # labels = ["PhEW,25/288", "PhEW,25/144"]

    modelnames = ["l25n288-phew-m5","l25n144-phew-m5-spl"] # mi
    labels = ["PhEW,25/288", "PhEW,25/288, Split"]

    # modelnames = ["l25n288-phewoff-fw","l25n144-phewoff"] # mi
    # labels = ["PhEWOff,25/288", "PhEWOff,25/144"]
    
    # titlestr = [r'$11.0 < \log(\frac{M_{vir}}{M_\odot}) < 11.5$',\
    #             r'$11.85 < \log(\frac{M_{vir}}{M_\odot}) < 12.15$',\
    #             r'$12.85 < \log(\frac{M_{vir}}{M_\odot}) < 13.15$']
    titlestr = ["low-mass",\
                "intermediate",\
                "massive"]
    lstyles = ["-", "--", ":"]
    moster18 = [8.92, 10.4, 10.96] # Moster 18 M* for the three bins
    moster18 = (array(moster18) - array([11.25, 12.0, 13.0])) - log10(OMEGAB/OMEGAM)
    from matplotlib.lines import Line2D
    lgds = []
    for fi in range(3): # Mass bins
        for mi in [0, 1]: # Models
            sfrinfoname, soname = find_fnames(modelnames[mi], mstrs[fi])
            print sfrinfoname, soname
            if(mi == -1):
                stars = load_central_stars(sfrinfoname, fformat="old")
                abins = build_abins(stars, cumulative=False)
            else:
                stars = load_central_stars(sfrinfoname, fformat="new")
                abins = build_abins_phew(stars, cumulative=False)            
            abins.totalmass = find_total_mass(stars, soname)
            print "Total Mass ----> ", abins.totalmass
            abins.draw(axs[3+fi], lstyle=lstyles[mi], method="differential") # Normed cumulative                
            abins.cumulate()
            abins.draw(axs[0+fi], lstyle=lstyles[mi], method="normed") # Normed by total Mass
            abins.draw(axs[6+fi], lstyle=lstyles[mi], method="cumulative") # Normed cumulative
            axs[fi].set_title(titlestr[fi])
            if(fi == 0):
                lgds.append(Line2D([0], [0], color="black", linestyle=lstyles[mi], label=labels[mi]))
        axs[fi].plot([log10(3.), log10(3.)], [-4.0, 0.0], ":", color="lightgrey")
        axs[fi].plot([log10(2.), log10(2.)], [-4.0, 0.0], ":", color="lightgrey")
        axs[fi].plot(log10(1.1), moster18[fi], "*", color="black", markersize=18)
        axs[fi+3].plot([log10(3.), log10(3.)], [-5.0, 0.0], ":", color="lightgrey")
        axs[fi+3].plot([log10(2.), log10(2.)], [-5.0, 0.0], ":", color="lightgrey")        
        axs[fi+6].plot([log10(3.), log10(3.)], [0.0, 1.1], ":", color="lightgrey")
        axs[fi+6].plot([log10(2.), log10(2.)], [0.0, 1.1], ":", color="lightgrey")        

    # from matplotlib.lines import Line2D
    # lgds = [Line2D([0], [0], color="black", linestyle=":", label="G3Cool"),\
    #         Line2D([0], [0], color="black", linestyle="-", label="GIZMO")]
    axs[0].legend(handles=lgds, loc="upper left")
    lgds = [Line2D([0], [0], color="black", linestyle="-", label="total"),\
            Line2D([0], [0], color="blue", linestyle="-", label="cold"),\
            Line2D([0], [0], color="red", linestyle="-", label="hot")]
    if(SEPARATE_COLD_HOT_WINDS):
        lgds.append(Line2D([0], [0], color=color_coldw, linestyle="-", label="cold wind"))
        lgds.append(Line2D([0], [0], color=color_hotw, linestyle="-", label="hot wind"))
    # lgds.append(Line2D([0], [0], color="green", linestyle="-", label="wind"))
    if(SHOW_WINDS_MIXED):
        lgds.append(Line2D([0], [0], color="cyan", linestyle="-", label="wmix"))
    lgds.append(Line2D([0], [0], color="lightgrey", linestyle="-", label="M(z)/M(z=0)"))
    axs[6].legend(handles=lgds, loc="upper left")

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

# fstarinfo, fso = find_fnames("l25n144-phew","mh12")
# stars = load_central_stars(fstarinfo)
