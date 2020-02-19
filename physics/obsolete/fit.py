# Try to find some relations between the Mach number and the acceleration.
# Not very successful

import ioformat
import matplotlib.pyplot as plt
from scipy import sqrt, array, log10, linspace
from astropy.table import Table
from astropy.io import ascii
from astroconst import pc,ac
import cosmology
import bowshock

import matplotlib as mpl
mpl.rcParams['mathtext.default'] = "default"
mpl.rcParams['axes.labelsize'] = "medium"

bs16 = Table.read("models_bs16.dat", format="ascii", data_start=0)
nmodels = len(bs16['ID'])
acc = array([0.0] * nmodels)
v0 = array([0.0] * nmodels)
p_ps = array([1.0] * nmodels) * pc.k / 300. * 3.e6
p_ps *= bowshock.rh_pressure_ratio(bs16['mach'])

fac = 1.375e12 # kTc/0.62m_h

clrlist = ['blue', 'green', 'orange', 'teal', 'purple', 'red', 'black']

for i in range(nmodels):
    acc[i] = (bs16['v25'][i]-bs16['v75'][i])/(bs16['t25'][i]-bs16['t75'][i])/bs16['tcc'][i] # km/s/Myr
    v0[i] = bs16['v50'][i] - acc[i] * bs16['tcc'][i] * bs16['t50'][i] # The initial velocity
    
for i in range(nmodels):
    plt.plot(bs16['mach'][i], acc[i], ".", color=clrlist[i], markersize=10)
    # plt.plot(p_ps, acc, ".")
    # plt.plot(bs16['mach'][i], v0[i], ".", color=clrlist[i], markersize=10)
plt.show()
