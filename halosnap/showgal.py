import matplotlib.pyplot as plt
from astroconst import pc, ac
from numpy import log10, logspace, inf
import gadget_units

fname = "/nas/astro-th/shuiyao/l25n144-test/gal_z108.stat"

GalTxt = "GIDX" # MSTAR, GID, GIDX
LABEL = False
USTC = False
XBOUNDS = (0., 25.)
YBOUNDS = (0., 25.)
LBOX = 25000. # kpc
NORM_SIZE = 0.5e11
HUBBLEPARAM = 0.7
XYRANGE = True
MGAL_MIN = 10.5
XRANGE = [0., 25.]
YRANGE = [0., 25.]

# units turns tipsy format output to comoving c.g.s units
units = gadget_units.tipsy_unit(LBOX, HUBBLEPARAM)

class parameters():
    def __init__(self):
        self.h0 = 0.70 # 0.72
        self.rho = 1.8791e-29 * self.h0 * self.h0
        self.l = LBOX # Mpc
        # For USTC runs: a 1000 factor is there
        if(USTC == True): # Note I added the self.h0 * self.h0 on 05/22/19
            # self.m = self.rho * (self.l * ac.mpc / 1000.) ** 3 / ac.msolar / (self.h0 * self.h0)
            self.m = self.rho * (self.l * ac.mpc) ** 3 / ac.msolar / (self.h0 * self.h0)            
        else:
            self.m = self.rho * (self.l * ac.mpc) ** 3 / ac.msolar / (self.h0 * self.h0)

class skid_group():
    def __init__(self):
        self.gid = 0
        self.mgal = 0.0
        self.mstar = 0.0
        self.mgas = 0.0
        self.pos = [0.0, 0.0, 0.0]
        self.vel = [0.0, 0.0, 0.0]        
        self.y = 0.0        
        self.size = 0.0

def read_gals(fname, units=units):
    f = open(fname, "r")
    gals = []
    for line in f:
        spt = line.split()
        gals.append(skid_group())
        gals[-1].gid = int(spt[0])
        gals[-1].mgal = float(spt[2]) * units.gm * 1.e10 # Solar mass
        gals[-1].mgas = float(spt[3]) * units.gm * 1.e10 # Solar mass
        gals[-1].mstar = float(spt[4]) * units.gm * 1.e10 # Solar mass
        gals[-1].pos[0] = (0.5 + float(spt[12])) * units.gl / 1.e3 # Mpc
        gals[-1].pos[1] = (0.5 + float(spt[13])) * units.gl / 1.e3 # Mpc
        gals[-1].pos[2] = (0.5 + float(spt[14])) * units.gl / 1.e3 # Mpc
        gals[-1].vel[0] = float(spt[15]) * units.gv # km/s
        gals[-1].vel[1] = float(spt[16]) * units.gv # km/s
        gals[-1].vel[2] = float(spt[17]) * units.gv # km/s       
        # gals[-1].size = 4. * log10(gals[-1].mstar) - 32.
        gals[-1].size = gals[-1].mgal / NORM_SIZE
        if(GalTxt == "MSTAR"):
            gals[-1].text = "%4.1f" % (log10(gals[-1].mstar))
        if(GalTxt == "GID"):        
            gals[-1].text = "%d" % (gals[-1].gid)
        if(GalTxt == "GIDX"):        
            gals[-1].text = "%d" % (len(gals)-1)
    f.close()
    return gals

def select_gals(gals, xlims=[], ylims=[], zlims=[], mlims=[]):
    newgals = []
    for gal in gals:
        if(len(xlims) == 2): # Upper and Lower limit set
            if(gal.pos[0] < xlims[0] or gal.pos[0] > xlims[1]): continue
        if(len(ylims) == 2): # Upper and Lower limit set
            if(gal.pos[1] < ylims[0] or gal.pos[1] > ylims[1]): continue
        if(len(zlims) == 2): # Upper and Lower limit set
            if(gal.pos[2] < zlims[0] or gal.pos[2] > zlims[1]): continue
        if(len(mlims) == 2): # Upper and Lower limit set
            if(gal.mstar < mlims[0] or gal.mstar > mlims[1]): continue
        newgals.append(gal)
    return newgals
            

def draw(gals, ax=None, pcolor="blue", marker="o", tshift=0.0):
    # Note here tshift must be in the correct unit that gal.pos = gal.vel * tshift
    # Normally gal.vel * t.gadget = [kpc] but gal.pos here is in [Mpc]
    print ("Start Drawing.")
    x, y, ms, mgal, sizes = [], [], [], [], []
    # xs, ys = [], []
    if(ax == None):
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111)
    n = 0
    for gal in gals:
        if(XYRANGE == True):
            if(gal.pos[0] < XRANGE[0] or gal.pos[0] > XRANGE[1]): continue
            if(gal.pos[1] < YRANGE[0] or gal.pos[1] > YRANGE[1]): continue
        if(log10(gal.mgal) > MGAL_MIN):
            x.append(gal.pos[0] + gal.vel[0] * tshift)
            y.append(gal.pos[1] + gal.vel[1] * tshift)
            sizes.append(gal.size)
            # if(gal.mstar == 0):
            #     xs.append(gal.pos[0])
            #     ys.append(gal.pos[1])
            n += 1
    print ("Number of Galaxies Shown: ", n)
    ax.scatter(x, y, s=sizes, alpha=0.6, color=pcolor, marker=marker)
    # ax.plot(xs, ys, "r.", markersize=5)
    ax.set_xlim(XBOUNDS)
    ax.set_ylim(YBOUNDS)
    ax.set_xlabel("Mpc/h")
    ax.set_ylabel("Mpc/h")    

gals = read_gals(fname)
draw(gals)
