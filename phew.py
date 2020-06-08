from numpy import genfromtxt
from numpy import array, sqrt, linspace, log, log10, isnan
import matplotlib as mpl
import matplotlib.pyplot as plt
from cosmology import tcosmic
from astroconst import pc, ac
from random import random
import config_mpl

GAMMA = 5./3.
MASS_CLOUD = 2.0e38
CMAP = "jet"
NSKIP = 1

modelname = "norecouple"
model = "l25n144-phew-"+modelname
filename = "/proj/shuiyao/"+model+"/WINDS/z1/sorted.phews"
fphewsname = "/scratch/shuiyao/scidata/newwind/"+model+"/phewsinfo.z1"    

class PhEWFields():
    '''
    Data structure that defines the format of PhEW tracking files.

    Example: 
    >>> Fields = PhEWFields()
    >>> info_fields = Fields.get_field_info(["dv","f_dis","Key"], verbose=True)
    >>> tab = read_phew_tracks("phews.sample", info_fields)
    '''
    def __init__(self, fname="fields.dat"):
        self.fields = genfromtxt(fname, dtype="i8,S2,S10", names=True)
        self.idx = {}
        self.nfields = len(self.fields)
        for i in range(self.nfields): self.idx[self.fields['name'][i]] = self.fields['idx'][i]
    def get_field_info(self, names, verbose=False):
        '''
        names: A list of names to read from the phews.* file
        '''
        cols = []
        namelst = []
        dtypestr = ""
        for name in names:
            idx = self.idx[name]
            cols.append(idx)
            namelst.append(name)
            dtypestr += self.fields['dtype'][idx] + ","
        return (cols, dtypestr, namelst)
        if(verbose == True):
            print "idx: ", cols
            print "names: ", namelst
            print "dtypes: ", dtypestr

class PhEWParticle():
    def __init__(self, key):
        self.key = key
        self.color = "black"
        self.mvir = -1.0
        self.msub = -1.0
        self.rvir = -1.0
        self.machi = -1.0
        self.machr = -1.0
        self.atimei = -1.0
        self.atimer = -1.0
        self.mcloud = -1.0
        self.track = []

def read_phew_fields(fname, fields):
    tab = genfromtxt(fname, usecols=fields[0], dtype=fields[1], names=fields[2], skip_header=1)
    return tab

print "Compiled."

def create_phew_particles(PhEWTracks):
    PhEWParticles = []
    thisKey = -1
    for i in range(len(PhEWTracks)):
        if(PhEWTracks[i]['Key'] == thisKey):
            PhEWP.track.append(PhEWTracks[i])
        elif(PhEWTracks[i]['Key'] > thisKey):
            if(thisKey != -1):
                PhEWP.track = array(PhEWP.track, dtype=PhEWTracks.dtype)
                PhEWParticles.append(PhEWP)
            thisKey = PhEWTracks[i]['Key']
            PhEWP = PhEWParticle(thisKey) # Create a new one
            PhEWP.track.append(PhEWTracks[i])            
        else:
            raise ValueError, "Tracks data not sorted?"
    PhEWP.track = array(PhEWP.track, dtype=PhEWTracks.dtype)
    PhEWParticles.append(PhEWP)
    print "---> Created %d PhEW particles from tracks.\n" % (len(PhEWParticles))
    return PhEWParticles

def set_colors_phew_particles(PhEWParticles, field, vmin, vmax, logscale=False, cmap=plt.get_cmap(CMAP)):
    logvmin = log(vmin)
    if(logscale == False):
        color_grad = vmax - vmin
    else:
        color_grad = log(vmax) - log(vmin)
    for PhEWP in PhEWParticles:
        if(field == "Mvir"):
            v = PhEWP.mvir
        else:
            v = PhEWP.track[field][0]

        if(logscale == False):
            color_val = (v - vmin) / color_grad
        else:
            color_val = (log(v) - logvmin) / color_grad
        if(v < vmin): color_val = 0.0
        if(v > vmax): color_val = 1.0
        PhEWP.color = cmap(color_val)

counter_spurious_particles = [0, 0, 0, 0, 0, 0]        
def remove_spurious_particles(track):
    dr = track['dr']
    atime = track['atime']
    vrel = track['vrel']
    mc = track['M_cloud']
    # if(max(dr) >= 1000.): return 0
    if(dr[0] >= 50.): return 1
    for i in range(len(dr)-1):
        if(atime[i+1] - atime[i] > 0.01): return 2
        if(abs(dr[i+1] - dr[i]) > 0.2 * dr[-1]): return 3
        if(abs(vrel[i+1] - vrel[i]) > 0.2 * vrel[0]): return 4
        if(mc[i+1] > mc[i]): return 5
    return 0

def get_initwinds_and_rejoin_info(PhEWParticles, filename, fformat='New'):
    # Note: some particles has HID = 0;
    # Their halo properties are from the last line of sovcirc. (Rvir < 0)
    if(fformat == "New"):
        tab = genfromtxt(filename, names=True, dtype=('f8,f8,f8,f8,f8,f8,f8,i8,f8,f8,f8,i8'))
    else:
        tab = genfromtxt(filename, names=True, dtype=('f8,f8,f8,f8,f8,f8,f8,f8,f8,i8')) 
    key_to_idx = dict()
    for i in range(len(tab)):
        if(tab[i]['Rvir'] > 0): key_to_idx[tab[i]['PhEWKey']] = i
    for PhEWP in PhEWParticles:
        if(PhEWP.key in key_to_idx):
            idx = key_to_idx[PhEWP.key]
            cs_a = sqrt(GAMMA * pc.k * 10.**tab[idx]['T_a'] / (0.60 * pc.mh))
            PhEWP.atimei = tab[idx]['a_i']
            PhEWP.atimer = tab[idx]['a_rejoin']
            PhEWP.machi = tab[idx]['Vinit'] * 1.e5 / cs_a
            PhEWP.machr = tab[idx]['Mach']
            PhEWP.msub = tab[idx]['LogMsub']
            PhEWP.mvir = tab[idx]['LogMvir']
            PhEWP.rvir = tab[idx]['Rvir']
            PhEWP.mcloud = tab[idx]['M_c']

def select_hard_particles(PhEWParticles):
    selected = []
    for i, PhEWP in enumerate(PhEWParticles):
        if(PhEWP.atimer == -1.0):
            selected.append(PhEWP)
    return selected

def select_particles(PhEWParticles, ntot=60, mmin=11.0, mmax=13.5):
    selected = []
    nbins = ntot / 2
    mbins = [0] * nbins
    dm = (mmax - mmin) / (float)(nbins)
    for i, PhEWP in enumerate(PhEWParticles):
        idx = (PhEWP.mvir - mmin) / dm
        if(isnan(idx) or idx < 0.0 or idx >= nbins): continue
        idx = (int)(idx)
        if(mbins[idx] >= 2): continue
        mbins[idx] += 1
        selected.append(PhEWP)
    return selected

def select_random_particles(PhEWParticles):
    selected = []
    random_norm = 0.5
    for i, PhEWP in enumerate(PhEWParticles):    
        random_cut = random_norm + 0.3 * (PhEWP.mvir - 12.0)
        if(random() > random_cut or PhEWP.mvir >= 15.0):
            continue
        selected.append(PhEWP)
    return selected

def draw_phew_particles(PhEWParticles, ax, fieldx, fieldy, fieldc, \
                        nskip=1, logxscale=False, logyscale=False, \
                        color_min=None, color_max=None, logcscale=False, alpha=0.2, \
                        showidx=False, allparts=False):
    if(color_min == None): color_min = min(fieldc)
    if(color_max == None): color_max = max(fieldc)
    set_colors_phew_particles(PhEWParticles, fieldc, color_min, color_max, logscale=logcscale)
    iskip = 0
    print "------ ", fieldy, "------"
    for i, PhEWP in enumerate(PhEWParticles):
    #     if(PhEWP.mvir > 13.0): print i
            
        if(iskip < nskip):
            iskip += 1
            continue

        if(allparts == False):
            flag = remove_spurious_particles(PhEWP.track)
            counter_spurious_particles[flag] += 1
            if(flag): continue

        if(fieldc == "Mvir"): color_val = PhEWP.mvir
        else: color_val = PhEWP.track[fieldc][0]
        if(color_min < color_val < color_max):
            if(fieldx == "time" or fieldx == "dt"):
                xarr = tcosmic(PhEWP.track['atime']) - tcosmic(PhEWP.track['atime'][0])
                xarr /= 1.e6 # Myr
            else:
                xarr = PhEWP.track[fieldx]
            if(fieldy == 'r/rvir'):
                if(PhEWP.rvir <= 0): continue
                yarr = PhEWP.track['dr'] / PhEWP.rvir
            else:
                yarr = PhEWP.track[fieldy]
            ynorm = 1
            if(fieldy == 'M_cloud'): ynorm = MASS_CLOUD
            if(logxscale == True): xarr = log10(xarr)
            if(logyscale == True): yarr = log10(yarr)            
            ax.plot(xarr, yarr/ynorm, ".-", color=PhEWP.color, alpha=alpha, markersize=2)
            ax.plot(xarr[0], yarr[0]/ynorm, "^", color=PhEWP.color, markersize=6)
            ax.plot(xarr[-1], yarr[-1]/ynorm, "*", color=PhEWP.color, markersize=8)
            if(showidx):
                ax.text(xarr[0], yarr[0], str(i), fontsize=6, color=PhEWP.color)
        iskip = 0

# filename = "phews.sample"
# filename = "/proj/shuiyao/m6n64beta5/WINDS/z2/sorted.phews"
def load_particles(filename=filename, fphewsname=fphewsname, fformat='New'):
    print "Loading PhEW particles from ", fphewsname
    Fields = PhEWFields()
    info_fields = Fields.get_field_info(["atime","Mass","Key","dr","dv","vrel","Mach","M_cloud","rho_a", "rho_c", "Mstar","T_c","T_a","t_dis","t_evap","t_khi"], verbose=True)
    tab = read_phew_fields(filename, info_fields)
    sub_Tc = tab[tab['T_c'] < 1.e4]
    pparts = create_phew_particles(tab)
    get_initwinds_and_rejoin_info(pparts, fphewsname, fformat=fformat)
    return pparts

def draw_field():
    fig = plt.figure(1, figsize=(8,6))
    ax = fig.add_subplot(111)
    # draw_phew_particles(pparts, ax, 'dr', 'M_cloud', 'vrel', nskip=40, color_min=100.e5, color_max=900.e5)
    # draw_phew_particles(pparts, ax, 'dr', 'vrel', 'vrel', nskip=30, color_min=400.e5, color_max=900.e5)
    # draw_phew_particles(pparts, ax, 'dr', 'rho_a', 'vrel', nskip=5, logyscale=True, color_min=400.e5, color_max=900.e5, alpha=0.4)
    # draw_phew_particles(pparts, ax, 'dr', 'T_c', 'vrel', nskip=50, logyscale=True, color_min=50.e5, color_max=200.e5, alpha=0.4)
    # draw_phew_particles(pparts, ax, 'time', 'r/rvir', 'Mvir', nskip=30, logyscale=False, color_min=11.0, color_max=13.5, alpha=0.2, showidx=True)
    # draw_phew_particles(pparts, ax, 'time', 'Mass', 'vrel', nskip=10, logyscale=False, color_min=200.e5, color_max=900.e5, alpha=0.2)
    # draw_phew_particles(pparts, ax, 'atime', 'Mach', 'vrel', nskip=1, logyscale=False, color_min=500.e5, color_max=900.e5, alpha=0.2) # REJOIN
    parts_to_show = select_particles(pparts)    
    draw_phew_particles(parts_to_show, ax, 'time', 't_dis', 'Mvir', nskip=1, logyscale=False, color_min=11.0, color_max=13.5, alpha=0.5)
    print counter_spurious_particles
    ax.set_xlim(0, 500)    
    ax.set_ylim(0, 1000)
    ax.set_xlabel("Time [Myr]")
    ax.set_ylabel(r"$t_\mathrm{dis} [Myr]$")    
    plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
    plt.show()

def figure():
    from pltastro import frame, draw
    import config_mpl
    pparts = load_particles(filename, fphewsname)
    
    frm = frame.multi(3, 2)
    params = frm.params
    params.figsize = (9, 10)
    params.left = 0.12
    params.right = 0.95
    params.bottom = 0.28
    params.top = 0.93    
    params.wspace = 0.30
    params.hspace = 0.05
    panels = frm.panels
    # panels.set_xticks([0., 0.2, 0.4, 0.6, 0.8, 1.0])
    panels.set_xlabels(r"$Time [Myr]$")

    panels.set_xlims(0., 800.)
    panels.ylims[0] = (0., 1.)        
    panels.ylims[1] = (0., 20.)
    panels.ylims[2] = (0., 1.)
    # panels.ylims[3] = (3.5, 5.)
    panels.ylims[3] = (4.0, 7.0)    
    panels.ylims[4] = (-28., -22.)
    panels.ylims[5] = (-30., -24.)
    panels.yticks[0] = [0.2, 0.4, 0.6, 0.8, 1.0]
    panels.yticks[1] = [5., 10., 15., 20.]
    panels.yticks[2] = [0.2, 0.4, 0.6, 0.8, 1.0]    
    # panels.yticks[3] = [4.0, 4.5, 5.0]
    panels.yticks[3] = [5.0, 6.0, 7.0]    
    panels.yticks[4] = [-28., -26., -24., -22.]
    panels.yticks[5] = [-30., -28., -26., -24.]        
    for i in range(6): panels.yticksON[i] = True
    panels.ylabels[0] = r"$R/R_\mathrm{vir}$"
    panels.ylabels[1] = r"$\mathcal{M}$"
    panels.ylabels[2] = r"$M_\mathrm{c}$"
    panels.ylabels[3] = r"$T_\mathrm{a}$"
    panels.ylabels[4] = r"$\rho_\mathrm{c}$"
    panels.ylabels[5] = r"$\rho_\mathrm{a}$"            
    fig, axs = draw(frm)

    # parts_to_show = select_particles(pparts)
    parts_to_show = select_hard_particles(pparts)    
    draw_phew_particles(parts_to_show, axs[0], 'time', 'r/rvir', 'Mvir', nskip=NSKIP, logyscale=False, color_min=11.0, color_max=13.5, alpha=0.4)
    print counter_spurious_particles
    draw_phew_particles(parts_to_show, axs[1], 'time', 'Mach', 'Mvir', nskip=NSKIP, logyscale=False, color_min=11.0, color_max=13.5, alpha=0.4)
    draw_phew_particles(parts_to_show, axs[2], 'time', 'M_cloud', 'Mvir', nskip=NSKIP, logyscale=False, color_min=11.0, color_max=13.5, alpha=0.4)
    draw_phew_particles(parts_to_show, axs[3], 'time', 'T_a', 'Mvir', nskip=NSKIP, logyscale=True, color_min=11.0, color_max=13.5, alpha=0.4)
    draw_phew_particles(parts_to_show, axs[4], 'time', 'rho_c', 'Mvir', nskip=NSKIP, logyscale=True, color_min=11.0, color_max=13.5, alpha=0.4)
    draw_phew_particles(parts_to_show, axs[5], 'time', 'rho_a', 'Mvir', nskip=NSKIP, logyscale=True, color_min=11.0, color_max=13.5, alpha=0.4)

    axs[1].text(0.5, 0.85, modelname, fontsize=12, transform=axs[1].transAxes)
    axcbar = fig.add_axes([0.15,0.15,0.7,0.01])
    norm1 = mpl.colors.Normalize(vmin=11.0, vmax=13.5)
    cdcbar = mpl.colorbar.ColorbarBase(axcbar, cmap=plt.get_cmap(CMAP), norm=norm1, orientation="horizontal")
    cdcbar.set_ticks([11, 11.5, 12, 12.5, 13.0, 13.5])
    cdcbar.set_ticklabels(["11","11.5","12","12.5","13","13.5"])
    cdcbar.set_label(r"$M_\mathrm{vir}$")
    plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
    plt.show()

# Check cmap
def check_cmap(cmapname="jet"):
    fig = plt.figure(1, figsize=(6,6))
    ax = fig.add_subplot(111)
    cmap=plt.get_cmap(cmapname)
    x = linspace(0.0, 1.0, 200)
    y = (x+2.0)**(1./3.)
    ax.scatter(x, y, color=cmap(x), s=6)
    plt.show()
    
    
