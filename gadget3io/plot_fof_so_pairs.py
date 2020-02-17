import matplotlib.pyplot as plt
from scipy import loadtxt, log10
import bin1d
import config_mpl

MDARK = 5.16e8
MRES = log10(MDARK * 64.0)
MLIM = log10(MDARK * 16.0)
ymin, ymax = 0.0, 1.5

dtype={'names': ('idx','nn','mfof','mvir','mstar','rmin','nngb'),\
       'formats': ('i4','i4','f4','f4','f4','f4','i2')}
pairs = loadtxt("fofpairs_z2.dat", dtype=dtype, skiprows=1)
# pairs = loadtxt("fofpairs.dat", dtype=dtype, skiprows=1)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7,6))
# plt.plot(pairs['mvir'], (pairs['mvir']-pairs['mfof']), 'b.', markersize=2)

from scipy import meshgrid, histogram2d, linspace
x, y = bin1d.subsample(pairs['mvir'], pairs['mvir']-pairs['mfof'], nonzero=True)
xbins = linspace(MLIM, 13.0, 40)
ybins = linspace(ymin, ymax, 40)
H, xedges, yedges = histogram2d(x, y, bins=[xbins, ybins]) # Need Linear Bin?
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
# cset = plt.contourf(H.T, 8, origin='lower', extent=extent, cmap=plt.cm.Purples)
cset = plt.contour(H.T, 8, origin='lower', extent=extent, cmap=plt.cm.autumn)

plt.axis([MLIM, 13.0, ymin, ymax])
x, y = bin1d.subsample(pairs['mvir'], pairs['mvir']-pairs['mfof'], nonzero=True)
s = bin1d.bin1d(x, y, nbins=30, logbins=False, bounds=0.68)
xline, yline = s.value, s.median
uline, lline = s.ubound, s.lbound
plt.plot(xline, yline, "b-")
plt.plot(xline, uline, "b--")
plt.plot(xline, lline, "b--")
plt.vlines(MRES, ymin, ymax, 'black', '--')
plt.plot([11.5, 12.5], [0.2, 0.2+1.0*0.15], 'b-')
plt.xlabel(r'$\log(M_{vir}/M_\odot)$')
plt.ylabel(r'$\log(M_{vir}/M_{fof})$')
plt.title('z = 2', fontsize=16)
ax.text(0.70, 0.95, "p50n288-fiducial", fontsize=16, transform=ax.transAxes)
plt.savefig("./figures/mfof.eps")
plt.show()

