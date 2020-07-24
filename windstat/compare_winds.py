import matplotlib.pyplot as plt
from numpy import genfromtxt, array, log10, linspace
import config_mpl
import ioformat

# f1 = "/proj/shuiyao/l25n144-phew-rcloud/WINDS/winds.45"
# f2 = "/proj/shuiyao/l25n288-phew-m5/WINDS.red/winds.45"

f1 = "/proj/shuiyao/P50N288s/p50n288fiducial/WINDS/winds.3"
f2 = "/proj/shuiyao/p50n576fi/WINDS/winds.3"

# tab1 = genfromtxt(f1, names=True)
# tab2 = genfromtxt(f2, names=True)

sig1 = ioformat.rcol(f1, [3], linestart=1)
sig2 = ioformat.rcol(f2, [3], linestart=1)

# step = 10
# plt.plot(tab1['a'][::step], tab1['Siggal'][::step], "b.", alpha=0.2)
# plt.plot(tab2['a'][::step], tab2['Siggal'][::step], "r.", alpha=0.2)

# plt.hist(tab1['Siggal'], bins=linspace(0., 300., 50), color="blue", histtype="step")
# plt.hist(tab2['Siggal'][::8], bins=linspace(0., 300., 50), color="red", histtype="step")

plt.hist(sig1, bins=linspace(0., 300., 50), color="blue", histtype="step")
plt.hist(sig2[::4], bins=linspace(0., 300., 50), color="red", histtype="step")
plt.ylabel("Wind Count")
plt.xlabel("sigma_gal [km/s]")
# plt.legend(["25/144", "25/288 (*1/8)"])
plt.legend(["G3,50/288", "G3,50/576 (*1/8)"])
plt.savefig("/scratch/shuiyao/figures/tmp.pdf")

# plt.hist(tab1['Mgal'], bins=linspace(8., 13.5, 50), color="blue", histtype="step")
# plt.hist(tab2['Mgal'][::8], bins=linspace(8., 13.5, 50), color="red", histtype="step")
plt.show()

