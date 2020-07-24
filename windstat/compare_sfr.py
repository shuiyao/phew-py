import matplotlib.pyplot as plt
from numpy import genfromtxt, array, log10, cumsum

sfr1 = "/proj/shuiyao/l25n144-phew-rcloud/sfr.txt"
sfr2 = "/proj/shuiyao/l25n288-phew-m5/sfr.txt"

tab1 = genfromtxt(sfr1, names=("a","total_sm","totsfrrate","rate_in_msunperyear","total_sum_mass_stars"))
tab2 = genfromtxt(sfr2, names=("a","total_sm","totsfrrate","rate_in_msunperyear","total_sum_mass_stars"))

fig, ax = plt.subplots(1,1)
#field = "total_sm"
field = "totsfrrate"
# field = "rate_in_msunperyear"
#field = "total_sum_mass_stars"
step = 1
XMIN, XMAX = 0.0, 1.0
plt.plot(tab1['a'][::step], tab1[field][::step], "b-")
plt.plot(tab2['a'][::step], tab2[field][::step], "r-")
# plt.plot(tab1['a'][::step], cumsum(tab1[field])[::step], "b-")
# plt.plot(tab2['a'][::step], cumsum(tab2[field])[::step], "r-")
ax.set_xlim(XMIN, XMAX)
# plt.yscale("log")
plt.show()
