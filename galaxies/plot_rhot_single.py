import matplotlib.pyplot as plt
from scipy import meshgrid, reshape, log, log10, pi, ndimage
import matplotlib as mpl
mpl.rcParams["mathtext.default"] = "tt"
mpl.rcParams["axes.labelsize"] = "large"

CONTLEVELS = 25

model = "l25n288-phew-m5-ngb64"
#model = "m25n256-simba"
#model = "l50n288-gizmo-phewoff"
gridfile = model+"/"+"mrhot_z0_108"
#gridfile = "rhot_c0003_128"

REDSHIFT = 0.0

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

#fig, axarr = plt.subplots(1, 2, sharey=True)
rhoth = log10(rho_thresh(REDSHIFT)+1.)
print "Density Thresh: ", rhoth
fig = plt.figure(1, figsize=(6,8))
fig.subplots_adjust(top=0.80)
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
plt.colorbar(ticks=[0,5,10], orientation="horizontal")
z2 = ndimage.gaussian_filter(z, sigma=1.0, order=0)
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
    plt.title("Phase Diagram: z = "+str(REDSHIFT), y=1.15)
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
plt.show()
