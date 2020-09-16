from mymod import *
from cosmology import acosmic, tcosmic
import h5py
import os

logz = linspace(0., 1., 201)
redz = 10. ** logz - 1.

PLOT_3D = True

# model = "l25n288-phew-m5-spl"
model = "l12n144-phew-movie"
haloid = 1927 # 12.0
haloid = 1142 # 12.3
ELEV, AZIM = -20, 15 # Viewing Angle
ELEV_EDGEON, AZIM_EDGEON = 90, 0 # Viewing Angle

RLIM = 0.8
RLIM_ZOOM = 0.15
halostr = ("00000"+str(haloid))[-5:]
snapi = 110
zstr = ("000"+str(snapi))[-3:]

fbase = "/scratch/shuiyao/scidata/halosnap/" + model + "/"
outname = fbase + "box/box_" + halostr + "_" + zstr
print "Reading: ", outname
fout = open(outname, "r")
hx, hy, hz, hmass, hrad = fout.readline().split()
hx = (float)(hx)
hy = (float)(hy)
hz = (float)(hz)
hrad = (float)(hrad)
fout.close()
gasp = genfromtxt(outname, names=True, skip_header=2)

if(PLOT_3D):
    fig = plt.figure(1, figsize=(9,6))
    from mpl_toolkits.mplot3d import Axes3D
    ax = fig.add_axes([0.0,0.0,0.665,1.0], projection="3d")
    ax.set_facecolor("black")
    ax._axis3don = False
    ax_edgeon = fig.add_axes([0.67,0.0,0.33,0.495], projection="3d")
    ax_edgeon.set_facecolor("black")
    ax_edgeon._axis3don = False    
    ax_faceon = fig.add_axes([0.67,0.505,0.33,0.495], projection="3d")
    ax_faceon.set_facecolor("black")
    ax_faceon._axis3don = False
    # ax = fig.add_subplot(111, projection="3d")
else:
    fig, ax = plt.subplots(1,1, figsize=(7,7))
    plt.gca().patch.set_facecolor("black")

cmap_logt=plt.get_cmap("plasma") # 4.0 - 7.0
cmap_phew=plt.get_cmap("Greens") # 0.0 - 1.0
# Get the color list
COLOR_STAR = (65./255., 105./255., 225./255., 1.0) # royal blue
clrs, sizes = [], []
phews = gasp[gasp['Mc'] > 0]
for gp in gasp:
    sizes.append(3)
    if (gp['SFflag'] > 0): clrs.append(COLOR_STAR)
    elif (gp['Mc'] > 0): # PhEW
        clrs.append(cmap_phew(gp['Mc']))
        sizes[-1] = 20. * (gp['Mc']) ** 2
    else:
        tval = (gp['logT'] - 4.0) / 2.5
        tval = max(0.0, tval)
        tval = min(1.0, tval)
        clrs.append(cmap_logt(tval))
    # sizes[-1] = gp['Hsml'] * gp['Hsml'] / 300.
    

#ax.scatter(phews['x'], phews['y'], marker='o', c=clrs, s=sizes)
if(PLOT_3D):
    # Main
    ax.scatter(gasp['x'], gasp['y'], gasp['z'], marker='o', c=clrs, s=sizes, edgecolors='none')
    txt = "z = %3.1f" % (redz[::-1][snapi])
    ax.text2D(0.08, 0.92, txt, fontsize=16, color="lightgrey", weight='heavy', transform=ax.transAxes)
    ax.set_xlim(hx - RLIM * hrad, hx + RLIM * hrad)
    ax.set_ylim(hy - RLIM * hrad, hy + RLIM * hrad)
    ax.set_zlim(hz - RLIM * hrad, hz + RLIM * hrad)
    for si in range(len(sizes)):
        if(sizes[si] != 3): sizes[si] *= 2.0
    # Face On
    ax_faceon.scatter(gasp['x'], gasp['y'], gasp['z'], marker='o', c=clrs, s=sizes, edgecolors='none')
    ax_faceon.set_xlim(hx - RLIM_ZOOM * hrad, hx + RLIM_ZOOM * hrad)
    ax_faceon.set_ylim(hy - RLIM_ZOOM * hrad, hy + RLIM_ZOOM * hrad)
    ax_faceon.set_zlim(hz - RLIM_ZOOM * hrad, hz + RLIM_ZOOM * hrad)
    ax_faceon.view_init(ELEV, AZIM)    
    # Edge On
    ax_edgeon.scatter(gasp['x'], gasp['y'], gasp['z'], marker='o', c=clrs, s=sizes, edgecolors='none')
    ax_edgeon.set_xlim(hx - RLIM_ZOOM * hrad, hx + RLIM_ZOOM * hrad)
    ax_edgeon.set_ylim(hy - RLIM_ZOOM * hrad, hy + RLIM_ZOOM * hrad)
    ax_edgeon.set_zlim(hz - RLIM_ZOOM * hrad, hz + RLIM_ZOOM * hrad)
    ax_edgeon.view_init(ELEV_EDGEON, AZIM_EDGEON)    
else:
    ax.scatter(gasp['x'], gasp['y'], marker='o', c=clrs, s=sizes)
    ax.set_xlim(hx - RLIM * hrad, hx + RLIM * hrad)
    ax.set_ylim(hy - RLIM * hrad, hy + RLIM * hrad)
    ax.view_init(ELEV, AZIM)
    
# plt.subplots_adjust(top=1.,bottom=0.,left=0.,right=1.)
plt.savefig("/scratch/shuiyao/figures/tmp.png")
plt.show()

print "DONE."

# from mpl_toolkits.mplot3d import Axes3D
# ax = fig.add_subplot(111, projection="3d")
# plt.gca().patch.set_facecolor("black")
# ax._axis3don = False
# ax.scatter(dxg, dyg, dzg, marker=".", alpha=0.4, c=cgaslist, cmap=plt.cm.cool, s=25, edgecolors="none", norm=(vmin, vmax))
# plt.subplots_adjust(top=1.,bottom=0.,left=0.,right=1.,hspace=0.,wspace=0.02)
# ax.plot(dxs, dys, dzs, "g*", markersize=6)
