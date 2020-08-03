from mymod import *
from cosmology import acosmic, tcosmic
import h5py
import os

logz = linspace(0., 1., 201)
redz = 10. ** logz - 1.

PLOT_3D = True

# model = "l25n288-phew-m5-spl"
model = "l12n144-phew-movie"
haloid = 1503 # 11.0
haloid = 2318 # 12.0
#haloid = 1972 # 12.0
#haloid = 251 # 13.1

RLIM = 0.5
halostr = ("00000"+str(haloid))[-5:]
snapi = 78
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
    fig = plt.figure(1, figsize=(7,7))
    from mpl_toolkits.mplot3d import Axes3D
    ax = fig.add_subplot(111, projection="3d")
    ax._axis3don = False
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
        sizes[-1] = 50. * (gp['Mc']) ** 2
    else:
        tval = (gp['logT'] - 4.0) / 2.5
        tval = max(0.0, tval)
        tval = min(1.0, tval)
        clrs.append(cmap_logt(tval))
    # sizes[-1] = gp['Hsml'] * gp['Hsml'] / 300.
    

#ax.scatter(phews['x'], phews['y'], marker='o', c=clrs, s=sizes)
if(PLOT_3D):
    ax.scatter(gasp['x'], gasp['y'], gasp['z'], marker='o', c=clrs, s=sizes, edgecolors='none')
    ax.set_zlim(hz - 1.0 * hrad, hz + 1.0 * hrad)
    txt = "z = %3.1f" % (redz[::-1][snapi])
    ax.text2D(0.08, 0.92, txt, fontsize=16, color="lightgrey", weight='heavy', transform=ax.transAxes)
else:
    ax.scatter(gasp['x'], gasp['y'], marker='o', c=clrs, s=sizes)

ax.set_xlim(hx - RLIM * hrad, hx + RLIM * hrad)
ax.set_ylim(hy - RLIM * hrad, hy + RLIM * hrad)

plt.subplots_adjust(top=1.,bottom=0.,left=0.,right=1.,hspace=0.,wspace=0.02)
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
