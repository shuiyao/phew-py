import matplotlib.pyplot as plt
from scipy import meshgrid, reshape, log, log10, pi, linspace, ndimage
import matplotlib as mpl
mpl.rcParams["mathtext.default"] = "tt"
mpl.rcParams["axes.labelsize"] = "large"

CONTLEVELS = linspace(-6.,1.,20)

#Z_solar = 0.0122
Z_solar = 0.0189

#model = "l25n144-phewoff"
model = "l25n144-phew-mach1"
gridfiles = [\
    model+"/"+"tabmet_z0_108"\
    ]

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
rhoth = log10(rho_thresh(0.)+1.)
print "Density Thresh: ", rhoth
fig = plt.figure(1, figsize=(6,8))
fig.subplots_adjust(top=0.80)
n_sub = 1
for fname in gridfiles:
    ax = plt.subplot(1,1,n_sub)
#    ax2 = ax.twiny()
    f = open(fname, "r")
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
#        z.append(log(float(f.readline().split()[0])+1.))
        spt = f.readline().split()
        m = float(spt[1])
        met = float(spt[2])
        if m > 0:
            met = met / m
        else: met = 0.
        if met == 0:
            z.append(-6)
        else:
            z.append(log10(met/Z_solar)) # METALS
    x, y = meshgrid(xnodes, ynodes)
    z = reshape(z, (ncells_x, ncells_y)).T
    z2 = ndimage.gaussian_filter(z, sigma=3., order=0)
    cont = plt.contourf(x, y, z2, CONTLEVELS, cmap=plt.get_cmap("Purples"))
    f.close()
    plt.colorbar(ticks=[-5,-3,-1, 1], orientation="horizontal")
    plt.contour(x, y, z2, cont.levels[4::2], colors="#393939", linestyles="solid")
    
#    ax2.callbacks.connect("xlim_changed", update_ax2)
    plt.axhline(5, xnodes[0], xnodes[-1], linestyle=":", color="black")
    plt.plot([rhoth, rhoth], [ynodes[0], ynodes[-1]], "k:")
       
    
#     plt.plot(xnodes, zhist, "k-")

    ax.set_xlim(xnodes[0], xnodes[-1])
    ax.set_ylim(ynodes[0], ynodes[-1])
#    plt.axis([xnodes[0], xnodes[-1], ynodes[0], ynodes[-1]])
    if (n_sub == 1):
        plt.ylabel(r'$Log(T)$')
        # plt.title("DESPH", y=1.15)
        plt.title("Zoom-in, z = 1", y=1.15)        
    elif(n_sub==2):
        plt.title("PESPH-GW-5xEsn", y=1.15)
    else:
        plt.title("PESPH-New", y=1.15)
    n_sub = n_sub + 1

    ax.set_xlabel(r'$Log(\rho/<\rho>)$')
#    plt.colorbar(orientation="horizontal")
    ax2 = ax.twiny()
    update_ax2(ax)
    ax2.set_xlabel(r'$Log(n_H) [cm^{-3}]$')


#axcbar = fig.add_axes([0.1,0.9,0.8,0.05])

#fig.subplots_adjust(vspace=0)
plt.show()
