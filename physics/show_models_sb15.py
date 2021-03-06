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
vmax = [120., 120., 250., 80., 210., 320.]
dvlbl = [20., 20., 50., 15., 40., 60.]

def t_life(mach, alpha=12.5):
    msq = mach * mach
    t_shock = 2.0 * mach / (msq + 1) ** (1./2.)
    t_exp = alpha * (msq + 1) ** (1./6.)
    print "Shock Time: %5.3f tcc" % (t_shock)
    print "Expansion Time: %5.3f tcc" % (t_exp)
    print "Lifetime: %5.3f tcc" % (t_shock + t_exp)
    return t_shock + t_exp

def t25fit(mach):
    return 6.0 * sqrt(1. + mach)

fig = plt.figure(1, figsize=(6,4))
ax = fig.add_subplot(111)
for i in range(6):
    # axs[idx].plot(ana['t75'][i], ana['v75'][i], color=clrs[ana['Group'][i]], marker=symlst[0], markersize=ms)
    mach = sb15['Mach'][i]
    t25pred = sb15['t25'][i]
    t25sb15 = t25fit(mach)
    t25 = t_life(mach, 11.5*0.75)
    ax.plot(mach, t25, marker=".", color=clrs[i], markersize=12)
    ax.plot(mach, t25pred, marker="*", color=clrs[i], markersize=12)
    ax.plot(mach, t25sb15, marker="+", color=clrs[i], markersize=12)    
ax.set_xlabel("Mach")
ax.set_ylabel("t25 [tcc]")
plt.show()


