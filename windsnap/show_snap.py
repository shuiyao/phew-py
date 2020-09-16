import matplotlib.pyplot as plt
import ioformat
from scipy import linspace, logspace, meshgrid, array, reshape, ndimage, log10, nonzero, sqrt, pi, ravel
from matplotlib.colors import LogNorm
from pylab import Circle, setp, Rectangle
from astroconst import pc, ac
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import gridspec
from cosmology import tcosmic

flag = 1 # PHEW
#flag = 0 # GW
HUBBLEPARAM = 0.7
ASCALE = 0.8 # z = 0.25; 0.833333 z = 0.2
RCCORR = 1.0 # SHOULD BE REMOVED THEN
#RCCORR = ASCALE * ASCALE # SHOULD BE REMOVED THEN
LSIZE = 1.0
# UNIT_L = 12000.0
# UNIT_M = 433697.735404 * (12./25.)**3
UNIT_L = 25000.0
UNIT_M = 433697.735404
UNIT_T = sqrt(8.0*pi/3.0) * ac.kpc * 1.e3 / (100.0 * HUBBLEPARAM * 1.e5)
UNIT_V = UNIT_L * ac.kpc / UNIT_T / sqrt(ASCALE)
NC = 1000. # Now use nc[i]
USE_RCLOUD_AS_SIZE = True
# RFAC = 10.0 # Boost covering area by a factor of RFAC * RFAC
RFAC = 2.3 # Boost covering area by a factor of RFAC * RFAC

print "UNIT_V = ", UNIT_V, "[km/s]"

modelname = "l25n144-phew-rcloud"
DESCRIPTION = "PhEW"
#modelname = "l25n144-phewoff"
#DESCRIPTION = "PhEWOff"
fbase = "/scratch/shuiyao/scidata/windsnap/" + modelname + "/"
fgrid = fbase + "halo_"+modelname+".grid"
fwind = fbase + "halo_"+modelname+".wind"

# Underlying density map
PLOT_CONTOUR = True
PLOT_PCOLOR = True
PLOT_IONS = True
COLORBAR_OVI = False

# winds:
PLOT_WINDS = True
PLOT_HSMOOTH = False
PLOT_LOS = False
XLOS, YLOS = 0.207, -0.0635
PRINT_ZW = False
CDMIN, CDMAX = 1.e10, 1.e19 # HI
CDMIN2, CDMAX2 = 11.0, 16.0
CONTLEVELS = logspace(CDMIN2, CDMAX2, 6) # OVI

# XMIN, XMAX = -0.5, 0.5
# YMIN, YMAX = -0.5, 0.5
# ZMIN, ZMAX = -0.045, -0.030
# CLMIN, CLMAX = 3.0, 7.0 # Temperature
CLMIN, CLMAX = 0.0, 5.0 # LastSFTime

XMIN, XMAX = 0.2, 0.4
YMIN, YMAX = 0.2, 0.4
# XMIN, XMAX = -0.1, 0.1
# YMIN, YMAX = 0.0, 0.2
# ZMIN, ZMAX = -1.0, 1.0
ZMIN, ZMAX = -0.15, 0.15

if(PLOT_IONS == False):
    m, mo = ioformat.rcol(fgrid, [0,1], linestart=2)
else:
    m, mo = ioformat.rcol(fgrid, [2,3], linestart=2)
    
f = open(fgrid, "r")
spt = f.readline().split()
spt = f.readline().split()
nx, ny = int(spt[0]), int(spt[1])
f.close()

sa = (UNIT_L * ASCALE / HUBBLEPARAM * ac.kpc) ** 2 / (float(nx) * float(ny))
fac_col = UNIT_M * 1.e10 / HUBBLEPARAM * ac.msolar / sa / pc.mh
m = array(m) * fac_col
mo = array(mo) * fac_col

if(PLOT_WINDS == True):
    if(flag == 0):
        xw, yw, zw, hw, rho, temp = ioformat.rcol(fwind, [0,1,2,5,6,7])
        wflag = []
        for x in xw: wflag.append(1)
    else:
        xw, yw, zw, hw, rc, temp, nc, mc, lastsftime, wflag = ioformat.rcol(fwind, [0,1,2,5,6,7,8,9,10,11])
        rc = array(rc) / UNIT_L / RCCORR # Originally, Rc is physical
    hw = array(hw) / UNIT_L
# xw = array(xw) * UNIT_L * ASCALE / HUBBLEPARAM
# yw = array(yw) * UNIT_L * ASCALE / HUBBLEPARAM
# zw = array(zw) * UNIT_L * ASCALE

xbins = linspace(-LSIZE/2., LSIZE/2., nx)
ybins = linspace(-LSIZE/2., LSIZE/2., ny)
xgrid, ygrid = meshgrid(xbins, ybins)
# m, mo = array(m), array(mo)
# m *= UNIT_M * 1.e10 * ac.msolar / (pc.mh * (LSIZE/float(nx)*ac.kpc)**2)
# mo *= UNIT_M * 1.e10 * ac.msolar / (pc.mh * (LSIZE/float(ny)*ac.kpc)**2)
m = m + min(m[nonzero(m)]);
m = reshape(m, (nx, ny)).T
mo = mo + min(mo[nonzero(mo)]);
mo = reshape(mo, (nx, ny)).T
m2 = ndimage.gaussian_filter(m, sigma=1.0, order=0)
mo2 = ndimage.gaussian_filter(mo, sigma=1.0, order=0)

fig = plt.figure(1, figsize=(8,8))
# gs = gridspec.GridSpec(2,1,height_ratios=[4,1])
# ax = plt.subplot(gs[0])
ax = fig.add_subplot(111)

#plt.pcolor(xgrid, ygrid, mo2, cmap="YlGnBu", norm=LogNorm(vmin=mo2.min(), vmax=mo2.max()))
if(PLOT_PCOLOR == True):
    plt.pcolor(xgrid, ygrid, m2, cmap="Purples", norm=LogNorm(vmin=CDMIN, vmax=CDMAX))
if(PLOT_IONS == False):
    print "Mass Range: ", m2.min(), m2.max()
    print "Oxygen Mass Range: ", mo2.min(), mo2.max()
else:
    print "MHI Range: ", m2.min(), m2.max()
    print "MOVI Range: ", mo2.min(), mo2.max()
#plt.contour(xgrid, ygrid, mo2, CONTLEVELS, colors="black")
if(PLOT_CONTOUR == True):
    plt.contour(xgrid, ygrid, mo2, CONTLEVELS, cmap=plt.get_cmap("copper"), norm=LogNorm())
# ax.scatter(xw, yw, s=hw, c="pink")
if(PLOT_WINDS == True):
    cmap=plt.get_cmap("jet")
    if(PLOT_HSMOOTH):
        for i in range(len(xw)):
            if(flag == 1):
                if(wflag[i] > 0 and (ZMIN < zw[i] < ZMAX)):
                    if(not((XMIN < xw[i] < XMAX) and (YMIN < yw[i] < YMAX))): continue
                    # ax.add_patch(Rectangle((xw[i]-hw[i], yw[i]-0.0001), 2.0*hw[i], 0.0002, fc="pink", ec="grey"))
                    ax.plot([xw[i]-hw[i], xw[i]+hw[i]], [yw[i], yw[i]], color="grey")

    count = 0
    for i in range(len(xw)): # The winds
        if(wflag[i] > 0 and (ZMIN < zw[i] < ZMAX)):
            if(not((XMIN < xw[i] < XMAX) and (YMIN < yw[i] < YMAX))): continue
            tout = (tcosmic(ASCALE) - tcosmic(abs(lastsftime[i]))) / 1.e9
            # clr = int((log10(temp[i]) - CLMIN) * 255. / (CLMAX - CLMIN))
            # clr = int(mc[i] * 255.)
            clr = int((tout - CLMIN) * 255. / (CLMAX - CLMIN))
            if(clr < 0): clr = 0
            if(clr > 255): clr = 255
            clr = cmap(clr)
            if(flag == 1):
                if(USE_RCLOUD_AS_SIZE):
                    # ax.add_patch(Circle((xw[i], yw[i]), radius=RFAC*rc[i]*sqrt(nc[i]), fc=clr, ec="k"))
                    ax.add_patch(Circle((xw[i], yw[i]), radius=RFAC*rc[i]*sqrt(nc[i]), color=clr, fill=False))
                else:
                    ax.add_patch(Circle((xw[i], yw[i]), radius=(XMAX-XMIN)/200.0, fc=clr, ec="k"))
            else:
                ax.add_patch(Circle((xw[i], yw[i]), radius=hw[i], fc=clr, ec="k"))                
            if(PRINT_ZW): plt.text(xw[i], yw[i], str(zw[i]), fontsize=8, color="black")
            count += 1
    print "Number of Winds Displayed: ", count
            
        # ax.add_patch(Circle((xw[i], yw[i]), radius=hw[i]/10.0, fc="pink", ec="k"))
    # ax.set_xlabel("x [kpc]")
    # ax.set_ylabel("y [kpc]")
    # titlestr = "ID="+num+", Mvir="+str(log10(Mvir))[:5]+", Rvir="+str(Rvir)[:5]
    if(PLOT_LOS):
        ax.plot(XLOS, YLOS, "+", markersize=12, color="black")
        ax.set_xlim(XLOS-0.008, XLOS+0.008)
        ax.set_ylim(YLOS-0.008, YLOS+0.008)

ax.set_xlim(XMIN, XMAX)
ax.set_ylim(YMIN, YMAX)
ax.set_xlabel("x")
ax.set_ylabel("y")
titlestr = "L="+str(UNIT_L)+" a="+str(ASCALE)+" "+DESCRIPTION
plt.title(titlestr)

fig.subplots_adjust(top=0.9, bottom=0.15, left=0.15, right=0.85)

# ColorBar: Cloud Temperature
# ---
# axcbar = fig.add_axes([0.9,0.1,0.01,0.8])                                     
# norm1 = mpl.colors.Normalize(vmin=CLMIN, vmax=CLMAX)
# cdcbar = mpl.colorbar.ColorbarBase(axcbar, cmap=plt.cm.jet, norm=norm1, orientation="vertical")
# cdcbar.set_ticks([4, 5, 6, 7])
# cdcbar.set_ticklabels(["4","5","6","7"])
# cdcbar.set_label("Tc")

# ColorBar: Winds
# ---
axcbar = fig.add_axes([0.9,0.1,0.01,0.8])
norm1 = mpl.colors.Normalize(vmin=CLMIN, vmax=CLMAX)
cdcbar = mpl.colorbar.ColorbarBase(axcbar, cmap=plt.cm.jet, norm=norm1, orientation="vertical")
cdcbar.set_ticks([0, 1, 2, 3, 4, 5])
cdcbar.set_ticklabels(["0", "1", "2", "3", "4", "5"])
cdcbar.set_label("time [Gyr]")

# ColorBar: OVI
# ---
if(PLOT_CONTOUR and COLORBAR_OVI):
    axcbar = fig.add_axes([0.9,0.1,0.01,0.8])
    norm1 = mpl.colors.Normalize(vmin=CDMIN2, vmax=CDMAX2)
    cdcbar = mpl.colorbar.ColorbarBase(axcbar, cmap=plt.cm.copper, norm=norm1, orientation="vertical")
    cdcbar.set_ticks([11, 12, 13, 14, 15, 16])
    cdcbar.set_ticklabels(["11", "12", "13", "14", "15", "16"])
    cdcbar.set_label("N_OVI")

# ColorBar: HI
# ---
axcbar = fig.add_axes([0.1,0.1,0.75,0.01])                                     
norm1 = mpl.colors.Normalize(vmin=log10(CDMIN), vmax=log10(CDMAX))
cdcbar = mpl.colorbar.ColorbarBase(axcbar, cmap=plt.cm.Purples, norm=norm1, orientation="horizontal")
cdcbar.set_ticks([11, 13, 15, 17, 19])
cdcbar.set_ticklabels(["11", "13", "15", "17", "19"])
cdcbar.set_label("N_HI")

# fig = plt.figure(2, figsize=(8,8))
# ax2 = fig.add_subplot(111, projection = "3d")
# ax2.scatter(xw, yw, array(zw) * UNIT_L, "o", s=rc*sqrt(1000.), color="blue")
# # for i in range(len(xw)):
# #     ax2.plot_surface(xw[i], yw[i], zw[i],  rstride=rc[i]*sqrt(1000.), color='b', linewidth=0, alpha=0.5)    
# ax2.plot([0.0],[0.0],[0.0], "k+", markersize=12)
# # ax2.set_xlim(-50.0, 50.0)
# # ax2.set_ylim(-50.0, 50.0)
# # ax2.set_zlim(-50.0, 50.0)
# ax2.set_xlabel("x")
# ax2.set_ylabel("y")
# ax2.set_zlabel("z")

plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
plt.show()
