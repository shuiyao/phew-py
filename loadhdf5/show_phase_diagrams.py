import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from numpy import genfromtxt, log, log10, sqrt, array, linspace
from scipy import meshgrid, reshape, pi, ndimage
import matplotlib as mpl
from cosmology import tcosmic

mpl.rcParams["mathtext.default"] = "tt"
mpl.rcParams["axes.labelsize"] = "large"
CONTLEVELS = 15
CONT_FLOOR = 1.e-4
REDSHIFT = 0.2
snapnum = 98
field = "LastSFTime"

# Search "Paint PhEW Particles"

fbase = "/nas/astro-th-nas/shuiyao/"
# modelname = "l50n288-phew-m5-spl"
snapstr = ("000"+str(snapnum))[-3:]
gridfile, fphewsname = [], []
models = ["l25n288-phew-m5-spl", "l50n288-phew-m5-spl"]
for model in models:
    gridfile.append("/home/shuiyao_umass_edu/sci/phew-py/data/"+model+"/"+"mrhot_phew_"+snapstr)
    folder = fbase + model + "/"
    fphewsname.append(folder + "snapshot_" + snapstr + ".phews")

if(field == "f_cloud"):
    vmin, vmax = 0.0, 1.0
    logscale = False
if(field == "LastSFTime"):
    vmin, vmax = 0, 5.0 # Gyr
    logscale = False
if(field == "f_wind"):
    vmin, vmax = 1.e-5, 0.1
    logscale = True

def rho_thresh(z, Om=0.3, OL=0.7):
    f_omega = Om * (1+z)**3
    f_omega = f_omega / (f_omega + OL)
    return 6.*pi*pi*(1. + 0.4093*(1./f_omega-1.)**0.9052) - 1.

def nh(x):
    rho_crit_baryon = 1.879e-29 * 0.7 * 0.7 * 0.044
    mh = 1.6733e-24
    return log10(10**x * rho_crit_baryon * 0.76 / mh)

def update_ax2(ax):
    ax2 = ax.twiny()
    x1, x2 = ax.get_xlim()
    ax2.set_xlim(nh(x1), nh(x2))
    ax2.figure.canvas.draw()
    ax2.set_xlabel(r'$Log(n_H) [cm^{-3}]$')

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
            v = (tcosmic(1./(1.+REDSHIFT)) - tcosmic(abs(v))) / 1.e9
        if(logscale == False):
            color_val = (v - vmin) / color_grad
        else:
            color_val = (log(v) - logvmin) / color_grad
        if(v < vmin): color_val = 0.0
        if(v > vmax): color_val = 1.0
        clrsarr.append(color_val)
    return clrsarr

def set_sizes_phew_particles(PhEWParticles):
    sizearr = []
    for PhEWP in PhEWParticles:
        v = abs(PhEWP['f_cloud'])
        size_val = v ** 1.5 * 50
        sizearr.append(size_val)
    return sizearr

def draw_model(ax, gridfile, fphewsname, field="LastSFTime"):
    f = open(gridfile, "r")
    xnodes, ynodes = [], []
    z = []
    zw = []
    spt = f.readline().split()
    ncells_x = int(spt[0])
    ncells_y = int(spt[1])
    ncells = ncells_x * ncells_y
    for i in range(ncells_x):
        xnodes.append(float(f.readline().split()[0]))
    for i in range(ncells_y):
        ynodes.append(float(f.readline().split()[0]))
    for i in range(ncells):
        spt = f.readline().split()
        # z.append(log(float(f.readline().split()[0])+1.))    
        z.append(log10(float(spt[1])+CONT_FLOOR))
        # zw.append(log10(float(spt[2])))
        # zw.append(float(spt[2]))
        zw.append(log10(float(spt[2])*float(spt[1])+CONT_FLOOR))
    x, y = meshgrid(xnodes, ynodes)
    z = reshape(z, (ncells_x, ncells_y)).T
    zw = reshape(zw, (ncells_x, ncells_y)).T

    z2 = ndimage.gaussian_filter(zw, sigma=1.5, order=0)
    # cont = plt.contourf(x, y, z2, CONTLEVELS, cmap=plt.get_cmap("Purples"))
    cont = ax.contourf(x, y, z2, CONTLEVELS, cmap=plt.get_cmap("Purples"))

    # ax.pcolor(xnodes, ynodes, zw, cmap=plt.cm.Purples, norm=LogNorm(vmin=zw.min(), vmax=zw.max()))
    # ax.pcolor(xnodes, ynodes, zw, cmap=plt.cm.Purples, vmin=0.0, vmax=1.0)

    # plt.colorbar(ticks=[0,5,10], orientation="horizontal")

    # z2 = ndimage.gaussian_filter(z, sigma=1.5, order=0)
    # cont2 = plt.contour(x, y, z2, 6, colors="black")

    f.close()
    #    ax2.callbacks.connect("xlim_changed", update_ax2)
    ax.axhline(5, xnodes[0], xnodes[-1], linestyle=":", color="black")
    ax.plot([rhoth, rhoth], [ynodes[0], ynodes[-1]], "k:")
    ax.set_xlim(xnodes[0], xnodes[-1])
    ax.set_ylim(ynodes[0], ynodes[-1])
    #    plt.axis([xnodes[0], xnodes[-1], ynodes[0], ynodes[-1]])
    ax.set_ylabel(r'$Log(T)$')
    # plt.title("z = "+str(REDSHIFT), y=1.15)
    ax.text(-1.,7.5, "WHIM", color="green")
    ax.text(5.,7.5, "Hot", color="red")
    ax.text(0.,3.2, "Diffuse", color="purple")
    ax.text(3.,3.2, "Condensed", color="teal")

    ax.set_xlabel(r'$Log(\rho/<\rho>)$')

    # ---- Paint PhEW Particles
    phews = genfromtxt(fphewsname, names=True)
    # phews_past = phews[phews['f_cloud'] < 0]
    # phews_cur = phews[phews['f_cloud'] > 0]
    # sphs = genfromtxt(fname_sph, names=True)

    # OBSOLETE: Particles that used to be PhEW
    # ----------------------------------------------------------------
    #step = 50
    # clrsarr = set_colors_phew_particles(phews_past, field, vmin, vmax, logscale=logscale, cmap=cmap)
    # sizearr = set_sizes_phew_particles(phews_past)
    # ax.scatter(phews_past['rhogcm3'][::step], phews_past['LogTK'][::step], marker="o", c=clrsarr[::step], s=sizearr[::step], alpha=0.6, cmap=cmap)
    # ----------------------------------------------------------------

    step = (int)(len(phews) / 400)
    phews_cur = phews[::step]
    clrsarr = set_colors_phew_particles(phews_cur, field, vmin, vmax, logscale=logscale, cmap=cmap)
    sizearr = set_sizes_phew_particles(phews_cur)
    ax.scatter(phews_cur['rhoa'], phews_cur['Ta'], marker="o", c=clrsarr, s=sizearr, alpha=0.6, cmap=cmap)
    ax.scatter(phews_cur['rhoc'], phews_cur['Tc'], marker="x", c=clrsarr, s=sizearr, alpha=0.6, cmap=cmap)

    # ax.text(0.6, 0.1, modelname, fontsize=12, transform=ax.transAxes)

    xarr = [6.1]*11
    yarr = linspace(6.0, 7.8, 11)
    s = linspace(0., 1., 11)
    s = s ** 1.5 * 50.0
    ax.scatter(xarr, yarr, s=s, color="black")
    for i in range(11)[1:]:
        txt = "%3.1f" % (0.1*i)
        ax.text(xarr[i]+0.3, yarr[i], txt, color="black", fontsize=8)

    update_ax2(ax)

fig, axs = plt.subplots(1,2, figsize=(9,6))
cmap = plt.get_cmap("jet")
rhoth = log10(rho_thresh(REDSHIFT)+1.)
print ("Density Thresh: ", rhoth)
fig.subplots_adjust(top=0.80, bottom=0.2, wspace=0.15, right=0.9)
for i in range(len(gridfile)):
    print ("Model: ", gridfile[i])
    draw_model(axs[i], gridfile[i], fphewsname[i])

logscale = False    
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

plt.savefig("/home/shuiyao_umass_edu/figures/tmp.pdf")
plt.show()
plt.close('all')
print ("DONE")
