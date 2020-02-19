import ioformat
import matplotlib.pyplot as plt
from scipy import sqrt, array, log10, linspace
from astropy.table import Table
from astropy.io import ascii
import cosmology
from matplotlib import gridspec
from pylab import setp
import bowshock

sb15 = Table.read("models_sb15.dat", format="ascii", data_start=0)

symlst = ["*", "^", ".", ".", ".", "."]
clrs = ["blue", "orange", "red", "green", "yellow", "purple", "brown"]
vmax = [50., 100., 195., 60., 150., 200.]
mach = [3.8, 6.5, 11.4, 3.5, 3.6, 11.4]

fig = plt.figure(1, figsize=(6,4))
ax = fig.add_subplot(111)
vs, dps, txts = [], [], []
for i in range(6):
    vpred = 15. * sqrt(mach[i] * mach[i] + 1.)
    vs.append(vpred)
    ax.plot(mach[i], vpred, "+", color=clrs[i], markersize=12)
    p, = ax.plot(mach[i], vmax[i], "*", color=clrs[i], markersize=12)
    dps.append(p)
    txts.append(sb15['Modelname'][i])
txts[-1] = "x300v3000b"
plt.legend(dps, txts, loc=4, fontsize=12)
ax.set_xlabel("Mach")
ax.set_ylabel("Vmax [km/s]")
ax.set_ylim(30., 250.)
ax.plot(mach, vs, "k:")
plt.show()


