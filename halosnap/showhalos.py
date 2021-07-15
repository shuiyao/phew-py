import matplotlib.pyplot as plt
from astroconst import pc, ac
from numpy import log10, logspace, inf, array
import ioformat
from matplotlib.patches import Circle

model = "l12n144-phew-movie"
zstr = "400"

galname = "/nas/astro-th/shuiyao/"+model+"/gal_z"+zstr+".stat"
soname = "/nas/astro-th/shuiyao/"+model+"/so_z"+zstr+".sovcirc"

GalTxt = "GIDX" # MSTAR, GID, GIDX
LABEL = True
USTC = False
LBOX = 12. # mpc
XBOUNDS = (0., LBOX)
YBOUNDS = (0., LBOX)
HUBBLEPARAM = 0.7
XYRANGE = True
MSUB_MIN = 11.0
XRANGE = [1., LBOX - 1.]
YRANGE = [1., LBOX - 1.]
ZRANGE = [1., LBOX - 1.]

def draw(ax=None, pcolor="blue", marker="o", tshift=0.0):
    x, y, z = ioformat.rcol(galname, [18, 19, 20])
    x = (array(x) + 0.5) * LBOX
    y = (array(y) + 0.5) * LBOX
    z = (array(z) + 0.5) * LBOX
    msub, rsub = ioformat.rcol(soname, [6, 7], linestart=1)

    print ("Start Drawing.")
    if(ax == None):
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111)
    n = 0
    # xs, ys, sizes = [], [], []
    for i in range(len(x)):
        if(msub[i] == 0): continue
        if(XYRANGE == True):
            if(x[i] < XRANGE[0] or x[i] > XRANGE[1]): continue
            if(y[i] < YRANGE[0] or y[i] > YRANGE[1]): continue
            if(z[i] < ZRANGE[0] or z[i] > ZRANGE[1]): continue            
        if(log10(msub[i] / HUBBLEPARAM) > MSUB_MIN):
            # xs.append(x[i])
            # ys.append(y[i])
            # sizes.append(rsub[i] / 1.e3) # Convert to kpc
            plt.gca().add_artist(Circle((x[i], y[i]), rsub[i]/1.e3, alpha=0.6, fc="blue"))
            if(LABEL == True):
                lbl = "%d %4.1f" % (i+1, log10(msub[i]/HUBBLEPARAM))
                plt.text(x[i], y[i], lbl, fontsize=6)
            # print xs[-1], ys[-1], sizes[-1]
            n += 1
    print ("Number of Galaxies Shown: ", n)
    # ax.scatter(x, y, s=sizes, alpha=0.6, color=pcolor, marker=marker)
    # ax.plot(xs, ys, "r.", markersize=5)
    ax.set_xlim(XBOUNDS)
    ax.set_ylim(YBOUNDS)
    ax.set_xlabel("Mpc/h")
    ax.set_ylabel("Mpc/h")    
    plt.savefig("/scratch/shuiyao/figures/tmp.png")

gals = read_gals(fname)
draw(gals)
print ("compiled.")
