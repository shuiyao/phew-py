import matplotlib.pyplot as plt
from numpy import genfromtxt, array, log10, linspace
import config_mpl
import ioformat

f1 = "/proj/shuiyao/l25n144-phew-rcloud/WINDS/winds.45"
f2 = "/proj/shuiyao/l25n288-phew-m5/WINDS.red/winds.45"

# f1 = "/proj/shuiyao/P50N288s/p50n288fiducial/WINDS/winds.3"
# f2 = "/proj/shuiyao/p50n576fi/WINDS/winds.3"

# tab1 = genfromtxt(f1, names=True)
# tab2 = genfromtxt(f2, names=True)

a1, sig1 = ioformat.rcol(f1, [0,3], linestart=1)
a2, sig2 = ioformat.rcol(f2, [0,3], linestart=1)

# step = 10
# plt.plot(tab1['a'][::step], tab1['Siggal'][::step], "b.", alpha=0.2)
# plt.plot(tab2['a'][::step], tab2['Siggal'][::step], "r.", alpha=0.2)

# plt.hist(tab1['Siggal'], bins=linspace(0., 300., 50), color="blue", histtype="step")
# plt.hist(tab2['Siggal'][::8], bins=linspace(0., 300., 50), color="red", histtype="step")

a1, a2 = array(a1), array(a2)
sig1, sig2 = array(sig1), array(sig2)
sig1b = sig1[a1 > 0.5]
sig2b = sig2[a2 > 0.5]

plt.hist(sig1, bins=linspace(0., 150., 50), color="blue", histtype="step")
plt.hist(sig2[::4], bins=linspace(0., 150., 50), color="red", histtype="step")
plt.hist(sig1b, bins=linspace(0., 150., 50), color="cyan", histtype="step")
plt.hist(sig2b[::4], bins=linspace(0., 150., 50), color="magenta", histtype="step")
plt.ylabel("Wind Count")
plt.xlabel("$\sigma_{gal}$ [km/s]")
plt.legend(["25/144", "25/288 (*1/8)","25/144, z > 1","25/288, z > 1"])
# plt.legend(["G3,50/288", "G3,50/576 (*1/8)"])
plt.savefig("/scratch/shuiyao/figures/tmp.pdf")

# plt.hist(tab1['Mgal'], bins=linspace(8., 13.5, 50), color="blue", histtype="step")
# plt.hist(tab2['Mgal'][::8], bins=linspace(8., 13.5, 50), color="red", histtype="step")
plt.show()

