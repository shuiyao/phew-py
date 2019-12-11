from numpy import genfromtxt
from numpy import array, sqrt, linspace, log, log10
import matplotlib as mpl
import matplotlib.pyplot as plt

class PhEWFields():
    '''
    Data structure that defines the format of PhEW tracking files.

    Example: 
    >>> Fields = PhEWFields()
    >>> info_fields = Fields.get_field_info(["dv","f_dis","ID"], verbose=True)
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
    def __init__(self, pid):
        self.pid = pid
        self.color = "black"
        self.track = []

def read_phew_tracks(fname, fields):
    tab = genfromtxt(fname, usecols=info_fields[0], dtype=info_fields[1], names=info_fields[2], skip_header=1)
    return tab

print "Compiled."

def create_phew_particles(PhEWTracks):
    PhEWParticles = []
    thisID = -1
    for i in range(len(PhEWTracks)):
        if(PhEWTracks[i]['ID'] == thisID):
            PhEWP.track.append(PhEWTracks[i])
        elif(PhEWTracks[i]['ID'] > thisID):
            if(thisID != -1):
                PhEWP.track = array(PhEWP.track, dtype=PhEWTracks.dtype)
                PhEWParticles.append(PhEWP)
            thisID = PhEWTracks[i]['ID']
            PhEWP = PhEWParticle(thisID) # Create a new one
            PhEWP.track.append(PhEWTracks[i])            
        else:
            raise ValueError, "Tracks data not sorted?"
    PhEWP.track = array(PhEWP.track, dtype=PhEWTracks.dtype)
    PhEWParticles.append(PhEWP)
    print "---> Created %d PhEW particles from tracks.\n" % (len(PhEWParticles))
    return PhEWParticles

def set_colors_phew_particles(PhEWParticles, field, vmin, vmax, logscale=False, cmap=plt.get_cmap("rainbow")):
    logvmin = log(vmin)
    if(logscale == False):
        color_grad = vmax - vmin
    else:
        color_grad = log(vmax) - log(vmin)
    for PhEWP in PhEWParticles:
        v = PhEWP.track[field][0]
        if(logscale == False):
            color_val = (v - vmin) / color_grad
        else:
            color_val = (log(v) - logvmin) / color_grad
        if(v < vmin): color_val = 0.0
        if(v > vmax): color_val = 1.0
        PhEWP.color = cmap(color_val)

def draw_phew_particles(PhEWParticles, ax, fieldx, fieldy, fieldc, \
                        nskip=1, logxscale=False, logyscale=False, \
                        color_min=None, color_max=None, logcscale=False, alpha=0.2):
    if(color_min == None): color_min = min(fieldc)
    if(color_max == None): color_max = max(fieldc)
    set_colors_phew_particles(PhEWParticles, fieldc, color_min, color_max, logscale=logcscale)
    iskip = 0
    for PhEWP in PhEWParticles:
        if(iskip < nskip):
            iskip += 1
            continue
        if(max(PhEWP.track['dr']) >= 1000.): continue
        if(PhEWP.track['dr'][0] >= 10.): continue
        if(color_min < PhEWP.track[fieldc][0] < color_max):
            xarr = PhEWP.track[fieldx]
            yarr = PhEWP.track[fieldy]
            if(logxscale == True): xarr = log10(xarr)
            if(logyscale == True): yarr = log10(yarr)            
            ax.plot(xarr, yarr, ".-", color=PhEWP.color, alpha=alpha)
            ax.plot(xarr[0], yarr[0], "^", color=PhEWP.color, markersize=8)
            ax.plot(xarr[-1], yarr[-1], "*", color=PhEWP.color, markersize=8)
        iskip = 0

# filename = "phews.sample"
# filename = "/proj/shuiyao/m6n64beta1/WINDS/z2/sorted.phews"
# Fields = PhEWFields()
#info_fields = Fields.get_field_info(["atime","ID","dr","dv","vrel","M_c","rho_a","Mstar","T_c","T_a"], verbose=True)
#tab = read_phew_tracks(filename, info_fields)
# sub_Tc = tab[tab['T_c'] < 1.e4]
#pparts = create_phew_particles(tab)        

fig = plt.figure(1, figsize=(6,6))
ax = fig.add_subplot(111)
# draw_phew_particles(pparts, ax, 'dr', 'M_c', 'vrel', nskip=40, color_min=100.e5, color_max=900.e5)
# draw_phew_particles(pparts, ax, 'dr', 'vrel', 'vrel', nskip=30, color_min=400.e5, color_max=900.e5)
# draw_phew_particles(pparts, ax, 'dr', 'rho_a', 'vrel', nskip=5, logyscale=True, color_min=400.e5, color_max=900.e5, alpha=0.4)
draw_phew_particles(pparts, ax, 'dr', 'T_a', 'vrel', nskip=50, logyscale=True, color_min=50.e5, color_max=200.e5, alpha=0.4)
plt.show()
# cmap=plt.get_cmap("rainbow")
# x = linspace(-1.0, 3.0, 200)
# y = (x+2.0)**(1./3.)
# ax.scatter(x, y, color=cmap(x))
    
    
