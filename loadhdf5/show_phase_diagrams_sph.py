import matplotlib.pyplot as plt
from numpy import genfromtxt, log, log10, sqrt, array, linspace
from scipy import meshgrid, reshape, pi, ndimage
import matplotlib as mpl
from cosmology import tcosmic

mpl.rcParams["mathtext.default"] = "tt"
mpl.rcParams["axes.labelsize"] = "large"
CONTLEVELS = 15
REDSHIFT = 0.2

# Search "Paint PhEW Particles"

fbase = "/proj/shuiyao/"
modelname = "m5kh30fs10"
model = "l25n144-phew-"+modelname
gridfile = "/scratch/shuiyao/sci/PHEW_TEST/"+model+"/"+"mrhot_z0_100"
folder = fbase + model + "/"
snapnum = 100
snapstr = ("000"+str(snapnum))[-3:]
fname = folder + "snapshot_" + snapstr + ".sph"

def rho_thresh(z, Om=0.3, OL=0.7):
    f_omega = Om * (1+z)**3
    f_omega = f_omega / (f_omega + OL)
    return 6.*pi*pi*(1. + 0.4093*(1./f_omega-1.)**0.9052) - 1.

def nh(x):
    rho_crit_baryon = 1.879e-29 * 0.7 * 0.7 * 0.044
    mh = 1.6733e-24
    return log10(10**x * rho_crit_baryon * 0.76 / mh)

def update_ax2(ax):
    x1, x2 = ax.get_xlim()
    ax2.set_xlim(nh(x1), nh(x2))
    ax2.figure.canvas.draw()

def set_colors_phew_particles(PhEWParticles, field, vmin, vmax, logscale=False, cmap=plt.get_cmap("jet")):
    clrsarr = []
    logvmin = log(vmin)
    if(logscale == False):
        color_grad = vmax - vmin
    else:
        color_grad = log(vmax) - log(vmin)
    for PhEWP in PhEWParticles:
        v = abs(PhEWP[field])
        if(field == "LastSFTime"):
            v = (tcosmic(0.5) - tcosmic(abs(v))) / 1.e9
        if(logscale == False):
            color_val = (v - vmin) / color_grad
        else:
            color_val = (log(v) - logvmin) / color_grad
        if(v < vmin): color_val = 0.0
        if(v > vmax): color_val = 1.0
        clrsarr.append(color_val)
    return clrsarr

fig, ax = plt.subplots(1,1, figsize=(6,8))
phews = genfromtxt(fname, names=True)
phews = phews[phews['f_wind']>0]
cmap = plt.get_cmap("jet")

rhoth = log10(rho_thresh(REDSHIFT)+1.)
print "Density Thresh: ", rhoth
fig.subplots_adjust(top=0.80, bottom=0.2)
n_sub = 1
ax = plt.subplot(1,1,n_sub)
f = open(gridfile, "r")
xnodes, ynodes = [], []
z = []
spt = f.readline().split()
ncells_x = int(spt[0])
ncells_y = int(spt[1])
ncells = ncells_x * ncells_y
for i in range(ncells_x):
    xnodes.append(float(f.readline().split()[0]))
for i in range(ncells_y):
    ynodes.append(float(f.readline().split()[0]))
for i in range(ncells):
#    z.append(int(f.readline().split()[0]))
    z.append(log(float(f.readline().split()[0])+1.))
x, y = meshgrid(xnodes, ynodes)
z = reshape(z, (ncells_x, ncells_y)).T
colormap = plt.get
cont = plt.contourf(x, y, z, CONTLEVELS, cmap=plt.get_cmap("Purples"))

# plt.colorbar(ticks=[0,5,10], orientation="horizontal")

z2 = ndimage.gaussian_filter(z, sigma=1.5, order=0)
cont2 = plt.contour(x, y, z2, 5, colors="black")
f.close()
#    ax2.callbacks.connect("xlim_changed", update_ax2)
plt.axhline(5, xnodes[0], xnodes[-1], linestyle=":", color="black")
plt.plot([rhoth, rhoth], [ynodes[0], ynodes[-1]], "k:")
ax.set_xlim(xnodes[0], xnodes[-1])
ax.set_ylim(ynodes[0], ynodes[-1])
#    plt.axis([xnodes[0], xnodes[-1], ynodes[0], ynodes[-1]])
if (n_sub == 1):
    plt.ylabel(r'$Log(T)$')
    plt.title(modelname+", z = "+str(REDSHIFT), y=1.15)
    ax.text(-1.,7.5, "WHIM", color="green")
    ax.text(5.,7.5, "Hot", color="red")
    ax.text(0.,3.2, "Diffuse", color="purple")
    ax.text(3.,3.2, "Condensed", color="teal")
elif(n_sub == 2):
    # plt.title("PESPH-HM12-AC", y=1.15)
    plt.title("PESPH-GW-NoCut", y=1.15)
else:
    plt.title("PESPH-GW-5xEsn", y=1.15)
n_sub = n_sub + 1

ax.set_xlabel(r'$Log(\rho/<\rho>)$')
ax2 = ax.twiny()
update_ax2(ax)
ax2.set_xlabel(r'$Log(n_H) [cm^{-3}]$')

# axcbar = fig.add_axes([0.1,0.9,0.8,0.05])
# fig.subplots_adjust(vspace=0)

# ---- Paint PhEW Particles

field = "f_wind"
vmin, vmax = 1.e-5, 0.1
logscale = True

step = 10
clrsarr = set_colors_phew_particles(phews, field, vmin, vmax, logscale=logscale, cmap=cmap)
ax.scatter(phews['rhogcm3'][::step], phews['LogTK'][::step], marker="o", c=clrsarr[::step], alpha=0.6, s=6, cmap=cmap)
# ax.text(0.6, 0.1, modelname, fontsize=12, transform=ax.transAxes())

axcbar = fig.add_axes([0.15,0.1,0.7,0.015])
if(logscale == False):
    norm1 = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
else:
    norm1 = mpl.colors.Normalize(vmin=log10(vmin), vmax=log10(vmax))    
cdcbar = mpl.colorbar.ColorbarBase(axcbar, cmap=cmap, norm=norm1, orientation="horizontal")
if(field == "f_wind"):
    cdcbar.set_ticks([-5., -4., -3., -2., -1.])
    cdcbar.set_ticklabels(["-5","-4","-3","-2","-1"])
    cdcbar.set_label(r"$\log(f_\mathrm{wind})$")
if(field == "LastSFTime"):
    cdcbar.set_ticks([1.0, 2.0, 3.0, 4.0, 5.0])
    cdcbar.set_ticklabels(["1.0","2.0","3.0","4.0","5.0"])
    cdcbar.set_label(r"$t_\mathrm{out}[\mathrm{Gyr}^{-1}]$")

# xarr = [6.1]*11
# yarr = linspace(6.0, 7.8, 11)
# s = array(range(10)) * 25. / 10.
# ax.scatter(xarr, yarr, s=s, color="black")

plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
plt.close('all')
print "DONE"
