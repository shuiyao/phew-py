# write_wind_features()
# write_mloss_info()
# V25Vc()
# figure()
# draw_phase_diagram_contour()

from numpy import genfromtxt
from numpy import array, sqrt, linspace, log, log10, median, logspace, meshgrid, isnan
from scipy import histogram2d, ndimage
from pylab import setp
import matplotlib as mpl
import matplotlib.pyplot as plt
import config_mpl
import matplotlib.patches as mpatches
from matplotlib.mlab import griddata
from matplotlib.colors import LogNorm
from cosmology import tcosmic
from astroconst import pc, ac
from random import random
import phew

# Borrowed from phew.py

GAMMA = 5./3.
CMAP = "jet"
NSKIP = 1

zi = 1

zs = [2.0, 1.0, 0.2]
zstrs = ["z2", "z1", "z0"]

REDSHIFT, zstr = zs[zi], zstrs[zi]

modelname = "m4"
model = "l25n288-phew-"+modelname
MC_INIT = 2.0e37
# filename = "/proj/shuiyao/"+model+"/WINDS/z1/sorted.tracks"
fphewsname = "/scratch/shuiyao/scidata/newwind/"+model+"/phewsinfo."+zstr
fwindsinfo = "/scratch/shuiyao/scidata/newwind/"+model+"/windsinfo."+zstr
fmlossinfo = "/scratch/shuiyao/scidata/newwind/"+model+"/mlossinfo."+zstr

gridfile = "/scratch/shuiyao/sci/PHEW_TEST/"+model+"/"+"mrhot_z0_100"

key_to_phew = dict();
filename_phews = "/proj/shuiyao/"+model+"/WINDS/"+zstr+"/sorted.phews"

class PhEWFields():
    '''
    Data structure that defines the format of PhEW tracking files.

    Example: 
    >>> Fields = PhEWFields()
    >>> info_fields = Fields.get_field_info(["dv","f_dis","Key"], verbose=True)
    >>> tab = read_phew_tracks("phews.sample", info_fields)
    '''
    def __init__(self, fname="fields_tracks.dat"):
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
        self.atimei = -1.0
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
        if(abs(dr[i+1] - dr[i]) > 0.1 * dr[-1]): return 3
        if(abs(vrel[i+1] - vrel[i]) > 0.1 * vrel[0]): return 4
        if(mc[i+1] > mc[i]): return 5
    return 0

def get_initwinds_and_rejoin_info(PhEWParticles, filename):
    # Note: some particles has HID = 0;
    # Their halo properties are from the last line of sovcirc. (Rvir < 0)
    tab = genfromtxt(filename, names=True, dtype=('f8,f8,f8,f8,f8,f8,i8,f8,f8,f8,i8'))
    key_to_idx = dict()
    for i in range(len(tab)):
        if(tab[i]['Rvir'] > 0): key_to_idx[tab[i]['PhEWKey']] = i
    for PhEWP in PhEWParticles:
        if(PhEWP.key in key_to_idx):
            idx = key_to_idx[PhEWP.key]
            cs_a = sqrt(GAMMA * pc.k * 10.**tab[idx]['T_a'] / (0.60 * pc.mh))
            PhEWP.atimei = tab[idx]['a_i']
            PhEWP.machi = tab[idx]['Vinit'] * 1.e5 / cs_a
            PhEWP.msub = tab[idx]['LogMsub']
            PhEWP.mvir = tab[idx]['LogMvir']
            PhEWP.rvir = tab[idx]['Rvir']

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
                        showidx=False, post_recouple=False, phew=True):
    if(color_min == None): color_min = min(fieldc)
    if(color_max == None): color_max = max(fieldc)
    set_colors_phew_particles(PhEWParticles, fieldc, color_min, color_max, logscale=logcscale)
    iskip = 0
    count = 0
    keys = []
    print "------ ", fieldy, "------"
    for i, PhEWP in enumerate(PhEWParticles):
    #     if(PhEWP.mvir > 13.0): print i
        # if(PhEWP.track['flag'][-1] != 0):
        #     continue
            
        if(iskip < nskip):
            iskip += 1
            continue

        # We don't remove spurious particles here
        # flag = remove_spurious_particles(PhEWP.track)
        # counter_spurious_particles[flag] += 1
        # if(flag): continue

        if(fieldc == "Mvir"): color_val = PhEWP.mvir
        else: color_val = PhEWP.track[fieldc][0]
        if(color_min < color_val < color_max):
            # keys.append(PhEWP.key)
            track = PhEWP.track[PhEWP.track['atime'] < 0.833333]
            count += 1
            if(fieldx == "time" or fieldx == "dt"):
                xarr = tcosmic(track['atime']) - tcosmic(track['atime'][0])
                xarr /= 1.e6 # Myr
            else:
                xarr = track[fieldx]
            if(fieldy == 'r/rvir'):
                if(PhEWP.rvir <= 0): continue
                yarr = track['dr'] / PhEWP.rvir
            else:
                yarr = track[fieldy]
            if(fieldy == 'T'):
                for j in range(len(yarr)):
                    if(yarr[j] > 10.0): yarr[j] -= 10.0
            if(fieldx == 'rho'):
                xarr = log10(track['rho'])

            # Find if there's a recoupling point
            if(phew):
                wflag = track['flag']
                i_recouple = -1
                for j in range(len(wflag)):
                    if(wflag[j] == 0):
                        i_recouple = j-1
                        break
            else:
                i_recouple = -1
            ynorm = 1
            if(fieldy == 'M_cloud'): ynorm = MASS_CLOUD
            if(logxscale == True): xarr = log10(xarr)
            if(logyscale == True): yarr = log10(yarr)
            if(post_recouple == False):
                ax.plot(xarr, yarr/ynorm, ".-", color=PhEWP.color, alpha=alpha, markersize=2)
                ax.plot(xarr[0], yarr[0]/ynorm, "^", color=PhEWP.color, markersize=6)
                if(i_recouple != -1): # recoupling point
                    ax.plot(xarr[-1], yarr[-1]/ynorm, "*", color=PhEWP.color, markersize=8)
                else: # annihilated
                    ax.plot(xarr[-1], yarr[-1]/ynorm, "x", color=PhEWP.color, markersize=6)
            else: # Only draw tracks after recouple
                if(i_recouple != -1):
                    # ax.plot([xarr[0],xarr[i_recouple]], [yarr[0]/ynorm,yarr[i_recouple]/ynorm], "--", color=PhEWP.color, alpha=alpha)
                    # ax.plot(xarr[0:i_recouple], yarr[0:i_recouple]/ynorm, ":", color=PhEWP.color, alpha=alpha)
                    ax.plot(xarr[i_recouple:], yarr[i_recouple:]/ynorm, ".-", color=PhEWP.color, alpha=alpha, markersize=2)
                    ax.plot(xarr[-1], yarr[-1]/ynorm, "*", color=PhEWP.color, markersize=8)
                    keys.append(PhEWP.key)
            if(i_recouple != -1):
                ax.plot(xarr[i_recouple], yarr[i_recouple]/ynorm, "o", color=PhEWP.color, markersize=8)                    
            if(showidx):
                ax.text(xarr[-1], yarr[-1], wflag[-1], fontsize=10, color=PhEWP.color)
            iskip = 0
    print "Number of PhEWs: ", count
    return keys

# filename = "phews.sample"
# filename = "/proj/shuiyao/m6n64beta5/WINDS/z2/sorted.phews"
def load_particles(phew=True):
    Fields = PhEWFields()
    if(phew == True):
        info_fields = Fields.get_field_info(["atime","Key","rho","T","flag"], verbose=True)
    else:
        info_fields = Fields.get_field_info(["atime","Key","dr","dv","Mgal","rho","T"], verbose=True)        
    tab = read_phew_fields(filename, info_fields)
    pparts = create_phew_particles(tab)
    if(phew == True):
        get_initwinds_and_rejoin_info(pparts, fphewsname)
    print "Load ", len(pparts), "PhEW Particles."
    return pparts

def draw_field(pparts, nskip=NSKIP, logyscale=False):
    fig = plt.figure(1, figsize=(6,6))
    ax = fig.add_subplot(111)
    # draw_phew_particles(pparts, ax, 'rho', 'T', 'Mvir', nskip=NSKIP, logyscale=False, color_min=11.0, color_max=13.5, alpha=0.3, showidx=True)
    # draw_phew_particles(pparts, ax, 'rho', 'T', 'Mgal', nskip=nskip, logyscale=False, color_min=11.0, color_max=13.5, alpha=0.3, showidx=False, phew=False)
    # draw_phew_particles(pparts, ax, 'dt', 'T', 'Mgal', nskip=nskip, logyscale=logyscale, color_min=11.0, color_max=13.5, alpha=0.3, showidx=False, phew=False)
    draw_phew_particles(pparts, ax, 'rho', 'T', 'dv', nskip=nskip, logyscale=logyscale, color_min=200.0, color_max=1000., alpha=0.3, showidx=False, phew=False)
    # print counter_spurious_particles
    plt.title("GIZMO")
    # plt.xlabel("Time [Myr]")
    plt.xlabel("log(Rho) [c.g.s]")    
    plt.ylabel("log(T)")
    plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
    plt.show()
    # return pparts

def draw_phase_diagram_contour(ax, gridfile=gridfile):
    from scipy import meshgrid, ndimage, reshape
    CONT_FLOOR = 1.e-4
    CONTLEVELS = 8
    unit_Density = 1.8791e-29 * 0.045;
    f = open(gridfile, "r")
    xnodes, ynodes = [], []
    z = []
    spt = f.readline().split()
    ncells_x = int(spt[0])
    ncells_y = int(spt[1])
    ncells = ncells_x * ncells_y
    for i in range(ncells_x):
        xnodes.append(float(f.readline().split()[0]))
        xnodes[-1] += log10(unit_Density)
    for i in range(ncells_y):
        ynodes.append(float(f.readline().split()[0]))
    for i in range(ncells):
        spt = f.readline().split()
        z.append(log10(float(spt[1])+CONT_FLOOR))
    x, y = meshgrid(xnodes, ynodes)
    z = reshape(z, (ncells_x, ncells_y)).T
    z2 = ndimage.gaussian_filter(z, sigma=2.0, order=0)
    # cont = ax.contour(x, y, z2, 5, colors="lightgrey", linestyles="dashed")
    cont = ax.contourf(x, y, z2, CONTLEVELS, cmap=plt.get_cmap("Purples")) 
    f.close()

def figure():
    from pltastro import frame, draw
    import config_mpl
    frm = frame.multi(2, 1)
    params = frm.params
    params.figsize = (5, 9)
    params.hspace = 0.35
    params.top = 0.93
    params.bottom = 0.28
    params.left = 0.20
    panels = frm.panels
    panels.set_ylims(2.5, 8.)
    panels.xlims[0] = (0., 6000.)
    panels.xlims[1] = (-31., -23.)
    panels.xticks[0] = [0, 2000, 4000, 6000]
    panels.xticks[1] = [-31., -29., -27., -25., -23.]
    panels.set_yticks([3, 4, 5, 6, 7])
    panels.ylabels[0] = r'$\log(T)$'
    panels.ylabels[1] = r'$\log(T)$'
    panels.xlabels[0] = "Time [Myr]"
    panels.xlabels[1] = r"$\log(\rho_a\ [\mathrm{g/cm^3}])$"
    for i in range(2): panels.yticksON[i] = True
    for i in range(2): panels.xticksON[i] = True

    fig, axs = draw(frm)
    
    pparts = load_particles() # tracks
    print "Draw Tracks....."
    keys = draw_phew_particles(pparts, axs[0], 'time', 'T', 'Mvir', nskip=NSKIP, logyscale=False, color_min=11.0, color_max=13.5, alpha=0.3, post_recouple=False)    
    keys = draw_phew_particles(pparts, axs[1], 'rho', 'T', 'Mvir', nskip=NSKIP, logyscale=False, color_min=11.0, color_max=13.5, alpha=0.3, post_recouple=False)

    phews = phew.load_particles(filename_phews, fphewsname) # phews
    for p in phews: key_to_phew[p.key] = p
    phews_to_show = []
    for key in keys: phews_to_show.append(key_to_phew[key])
    print "Draw PhEWs (%d) ....." % (len(phews_to_show))
    phew.draw_phew_particles(phews_to_show, axs[0], 'time', 'T_c', 'Mvir', nskip=1, logyscale=True, logxscale=False, color_min=11.0, color_max=13.5, alpha=0.3, allparts=True)
    phew.draw_phew_particles(phews_to_show, axs[1], 'rho_a', 'T_c', 'Mvir', nskip=1, logyscale=True, logxscale=True, color_min=11.0, color_max=13.5, alpha=0.3, allparts=True)

    draw_phase_diagram_contour(axs[1])

    # axs[0].text(0.5, 0.85, modelname, fontsize=12, transform=axs[1].transAxes)
    axs[0].set_title(model+", z = "+str(REDSHIFT))
    
    axcbar = fig.add_axes([0.15,0.15,0.7,0.01])
    norm1 = mpl.colors.Normalize(vmin=11.0, vmax=13.5)
    cdcbar = mpl.colorbar.ColorbarBase(axcbar, cmap=plt.get_cmap(CMAP), norm=norm1, orientation="horizontal")
    cdcbar.set_ticks([11, 11.5, 12, 12.5, 13.0, 13.5])
    cdcbar.set_ticklabels(["11","11.5","12","12.5","13","13.5"])
    cdcbar.set_label(r"$M_\mathrm{vir}$")
    
    plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
    plt.show()


def write_wind_features(foutname=fwindsinfo):
    PhEWParticles = phew.load_particles(filename_phews, fphewsname) # phews    
    fout = open(foutname, "w")
    fout.write("#ID Mvir Rvir Vinit V25 Rreturn\n")
    count, writecount = 0, 0
    for n, PhEWP in enumerate(PhEWParticles):
        flag = remove_spurious_particles(PhEWP.track)
        counter_spurious_particles[flag] += 1
        if(flag == 5): print n
        if(flag): continue
        PhEWP.mvir = 10. ** PhEWP.mvir
        vc, v25, vinit = 0., 0., 0.            
        x, y = [PhEWP.track['dr'][0]], [PhEWP.track['dv'][0]]
        dtime = (tcosmic(PhEWP.track['atime']) - tcosmic(PhEWP.track['atime'][0])) / 1.e6  
        ddr, ddt = [], []
        i25, r25 = 0, -1
        vel25 = -1
        ireturn = 0
        count += 1
        vc = sqrt(pc.G * PhEWP.mvir * ac.msolar / (PhEWP.rvir * ac.kpc)) / 1.e5
        for i in range(len(PhEWP.track['dr']))[1:]:
            x.append(abs(PhEWP.track['dr'][i]))
            y.append(PhEWP.track['dv'][i])
            ddt.append(dtime[i] - dtime[i-1])
            ddr.append(x[i]-x[i-1])
            if(i25 == 0 and abs(PhEWP.track['dr'][i]) > PhEWP.rvir * 0.25 and dtime[i]<1000.): # WARNING: dt < 400.
                i25 = i
            if(ireturn == 0 and ((ddr[i-1] < 0. and PhEWP.track['dv'][i] < vc) or ddt[i-1]>30.)):
                ireturn = i
                rreturn = x[i]
        if(ireturn == 0):
            ireturn = -1
            rreturn = -1
        if(len(ddr) == 0):
            continue
        if(PhEWP.rvir == -1 or isnan(PhEWP.mvir)):
            print "Warning: ", PhEWP.key, PhEWP.mvir, PhEWP.rvir
            continue
        if(max(ddr) < 100. and min(ddr) > -100.):
            writecount += 1
            if(i25 > 0):
                v25 = y[i25]
            else:
                v25 = -1
            vinit = y[0]
            # fout.write("#ID, Mvir, Vc, Vinit, V25\n")            
            outstr = str(PhEWP.key)+" "+str(log10(PhEWP.mvir))+" "+str(PhEWP.rvir)
            outstr += " "+str(vinit)+" "+str(v25)+" "+str(rreturn)
            outstr += "\n"
            fout.write(outstr)
    print counter_spurious_particles            
    print "Written, Count, All: ", writecount, count, len(PhEWParticles)
    fout.close()

def write_mloss_info(foutname=fmlossinfo):
    PhEWParticles = phew.load_particles(filename_phews, fphewsname) # phews    
    fout = open(foutname, "w")
    fout.write("#ID Mvir Rvir R75 R50 R25 Rlast t75 t50 t25 tlast Mach75 Mach50 Mach25 Machlast\n")
    for PhEWP in PhEWParticles:
        flag = remove_spurious_particles(PhEWP.track)
        counter_spurious_particles[flag] += 1
        # if(flag): continue
        if(isnan(PhEWP.mvir) or PhEWP.mvir < 0): continue
        r25, r50, r75, rlast = -1., -1., -1., -1.
        mach25, mach50, mach75, machlast = -1., -1., -1., -1.
        t25, t50, t75, tlast = -1., -1., -1., -1.        
        PhEWP.track['M_cloud'] /= MC_INIT
        dtime = (tcosmic(PhEWP.track['atime']) - tcosmic(PhEWP.track['atime'][0])) / 1.e6          
        for i in range(len(PhEWP.track['dr']))[1:]:
            if(r75 < 0.0 and PhEWP.track['M_cloud'][i] <= 0.75):
                r75, mach75, t75 = abs(PhEWP.track['dr'][i]), PhEWP.track['Mach'][i], dtime[i]
            if(r50 < 0.0 and PhEWP.track['M_cloud'][i] <= 0.50):
                r50, mach50, t50 = abs(PhEWP.track['dr'][i]), PhEWP.track['Mach'][i], dtime[i]
            if(r25 < 0.0 and PhEWP.track['M_cloud'][i] <= 0.25):
                r25, mach25, t25 = abs(PhEWP.track['dr'][i]), PhEWP.track['Mach'][i], dtime[i]
        if(r25 > 0):
            rlast, machlast, tlast = abs(PhEWP.track['dr'][-1]), PhEWP.track['Mach'][-1], dtime[i]
        outstr = "%10d %6.3f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f %6.3f %6.3f %6.3f\n" % \
                 (PhEWP.key, PhEWP.mvir, PhEWP.rvir, \
                  r75, r50, r25, rlast, \
                  t75, t50, t25, tlast,
                  mach75, mach50, mach25, machlast)
        fout.write(outstr)
    fout.close()

XMIN, XMAX = 30., 600.
YMIN, YMAX = 50., 2000.
XBINS, YBINS = 30., 30.
CONTLEVELS = 5
    
def plotmedian(x, y, clr, opt=1, nbins=10, ax=[], xmin=30., xmax=150):
    xmid, ymid, ylist = [], [], []
    dx = (xmax - xmin) / nbins
    for i in range(nbins):
        xmid.append(xmin + (0.5+i)*dx)
        ylist.append([])
        ymid.append(0.)
    for i in range(len(x)):
        idx = int((x[i] - xmin) / dx)
        if(idx < 0):
            idx = 0
        if(idx > nbins-1):
            idx = nbins-1
        ylist[idx].append(y[i])
    for i in range(nbins):
        if((len(ylist[i])==0 or max(ylist[i])<0.1) and i>0):
            ymid[i] = ymid[i-1]
            xmid[i] = xmid[i-1]
        else:
            ymid[i] = median(ylist[i])
    if(opt==1):
        ax.plot(xmid, ymid, ".-", color=clr)
    if(opt==2):
        plt.plot(xmid, ymid, ".-", color=clr)
    print xmid, ymid

def V25Vc(fname=fwindsinfo):
    import ioformat
    from matplotlib import gridspec
    
    fig = plt.figure(1, figsize=(8,8))
    gs = gridspec.GridSpec(2,1,height_ratios=[4,1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    mvir, rvir, vinit, v25, rreturn = ioformat.rcol(fname, [1,2,3,4,5], linestart=1)    
    vc = sqrt(pc.G * 10.**array(mvir) * ac.msolar / (array(rvir) * ac.kpc)) / 1.e5
    print "Mvir Range: ", min(mvir), max(mvir)
    print "Vc Range: ", min(vc), max(vc)
    # ax1.plot(vc[::100], v25[::100], "b.")
    # ax1.plot(vc2, v252, ".", color="teal")
    # Rreturn < 0.25 Rvir!
    for i in range(len(v25)):
        if(0.0 < rreturn[i] < 0.25*rvir[i]):
            v25[i] = -1
    plotmedian(vc, v25, "blue",nbins=15, ax=ax1, xmin=XMIN, xmax=XMAX)
    plotmedian(vc, vinit, "red",nbins=15, ax=ax1, xmin=XMIN, xmax=XMAX)
    xbins = logspace(log10(XMIN), log10(XMAX), XBINS)
    ybins = logspace(log10(YMIN), log10(YMAX), YBINS)
    xgrid, ygrid = meshgrid(xbins, ybins)
    z, edx, edy = histogram2d(vc, v25, bins=[xbins,ybins])
    z = z + 0.1
    z = z.T
    zf = ndimage.gaussian_filter(z, sigma=0.3, order=0)
    # cont = ax1.contour(xbins[1:], ybins[1:], zf, 10, colors="black")
    # cont = ax1.contourf(xbins[1:], ybins[1:], zf, CONTLEVELS, cmap=plt.cm.Purples, norm=LogNorm(vmin=z.min(), vmax=z.max()))
    # cont = ax1.contour(xbins[1:], ybins[1:], zf, CONTLEVELS, cmap=plt.cm.Reds, norm=LogNorm(vmin=z.min(), vmax=z.max()))    
    # cont = ax1.contourf(xbins[1:], ybins[1:], zf, 6, cmap=plt.cm.Purples, vmin=z.min(), vmax=z.max())
    # ax1.plot(vc[::200], v25[::200], "k.", markersize=2)
    ax1.pcolor(xbins[1:], ybins[1:], zf, cmap=plt.cm.Purples, norm=LogNorm(vmin=z.min(), vmax=z.max()))
    xline = linspace(XMIN, XMAX, 100)
    y50line = 0.854 * xline ** 1.12
    y95line = 1.85 * xline ** 1.10
    ax1.plot(xline, y50line, "k-")
    ax1.plot(xline, y95line, "k--")
    ax2.set_xlabel(r"$V_c [km/s]$")
    ax1.set_ylabel(r"$V_{25} [km/s]$")
    ax1.set_xlim(XMIN, XMAX)
    ax1.set_ylim(YMIN, YMAX)
    # ax1.set_ylim(9.,12.)
    ax1.set_title(modelname+", Z = "+str(REDSHIFT))
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.legend([r"$V_{25}$", r"$V_{init}$"], fontsize=16)
    fig.subplots_adjust(hspace=0)
    setp(ax1.get_xticklabels(),visible=False)
    # Count the fraction that makes NOT to R25
    FRAC_STEP = 25.
    frac = [0.] * int(XMAX/FRAC_STEP)
    Ncount = [0.] * int(XMAX/FRAC_STEP)
    x = linspace(10., XMAX, int(XMAX/FRAC_STEP)) - 5.
    for i in range(len(vc)):
        idx = int(vc[i]/FRAC_STEP)
        if(idx > int(XMAX/FRAC_STEP)-1):
            idx = int(XMAX/FRAC_STEP)-1
        Ncount[idx] += 1.
        if(v25[i] != -1):
            frac[idx] += 1.
        # if(FFORMAT == "NEW"):
        #     if(0. < rreturn[i] < 0.25 * rvir[i]):
        #         frac[idx] -= 1. # Rreturn < 0.25 Rvir!
    Ntot = sum(Ncount)
    for i in range(len(frac)):
        if(Ncount[i] > 0.):
            frac[i] = frac[i] / Ncount[i]
    ax2.plot(x, frac, "k.-")
    # ax2.plot(x, array(Ncount)/Ntot, "-", color="teal")
    Ncount_Norm = 2.*max(Ncount)/Ntot
    ax2.bar(x, (array(Ncount)/Ntot)/Ncount_Norm, align="center", width=0.8*(x[1]-x[0]), color="grey")
    ax2.set_xlim(XMIN,XMAX)
    ax2.set_ylim(0.,1.1)
    ax2.set_xscale("log")
    ax2.xaxis.set_ticks([XMIN, 50., 100., 200., XMAX])
    ax2.xaxis.set_ticklabels([str(XMIN), "50", "100", "200", str(XMAX)])
    ax2.yaxis.set_ticks([0., 0.2, 0.4, 0.6, 0.8])
    # ax2.set_ylabel(r"$f(R_{25} < R_{return})$")
    ax2.set_ylabel(r"$f(R_{25} < R_{recouple})$")    
    plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
    plt.show()

def scatterplot_vinit_Tfinal():    
    from pltastro import frame, draw, legend
    from numpy import histogram
    frm = frame.multi(2,1)
    pars = frm.params
    pars.left = 0.2
    pars.figsize = (6, 8)
    pnls = frm.panels
    pars.height_ratios = [4,1]
    pnls.ylabels = ["Vinit [km/s]", "frac"]
    pnls.set_xlabels("log(T_final)")
    pnls.set_xlims(3.0, 7.0)
    pnls.ylims[0] = (0.0, 1500.)
    pnls.ylims[1] = (0.0, 0.3)

    fig, axs = draw(frm)
    vi, Tf = [], []
    for p in pg3:
        vi.append(p.track['dv'][0])
        Tf.append(log10(p.track['T'][-1]))
        # Tf.append(p.track['T'][-1])    
    axs[0].plot(Tf, vi, "b.", alpha=0.3, markersize=4)
    hist, edges = histogram(Tf, bins=linspace(3, 7, 30))
    mid = 0.5*(edges[1:] + edges[:-1])
    norm = sum(hist)
    hist = array(hist) / (float)(norm)
    axs[1].plot(mid, hist, "b.-")
    vi, Tf = [], []
    for p in pparts:
        vi.append(p.track['dv'][0])
        Tf.append(p.track['T'][-1])
    axs[0].plot(Tf, vi, "r.", alpha=0.3, markersize=4)
    hist, edges = histogram(Tf, bins=linspace(3, 7, 30))
    mid = 0.5*(edges[1:] + edges[:-1])
    norm = sum(hist)
    hist = array(hist) / (float)(norm)
    axs[1].plot(mid, hist, "r.-")

    lgd = legend.legend(axs[0])
    lgd.addLine(("Gadget3", "blue", "-", 1))
    lgd.addLine(("GIZMO", "red", "-", 1))
    lgd.loc = "upper left"
    lgd.fontsize = 12
    lgd.draw()

    plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
    plt.show()
