from numpy import genfromtxt, logspace, linspace, median, log10, histogram, cumsum
from astroconst import pc, ac
import matplotlib.pyplot as plt
import config_mpl
from pltastro import legend
import pltastro

Mgas_lres = 0.00653561
Mgas_hres = 0.00081695

def reverse_cumsum(hist, edges):
    cen = 0.5 * (edges[1:] + edges[:-1])
    cumfrac = cumsum(hist[::-1])[::-1]
    cumfrac = cumfrac / (float)(cumfrac[0])
    return cen, cumfrac

SHOW_WIND_FRAC = False

# ---------------- STARINFO ----------------

# folder = "/scratch/shuiyao/scidata/gadget3io/l25n288-phew-m5/"
# fname = folder + "l25n288-phew-m5_108.starinfo"
# # folder = "/scratch/shuiyao/scidata/gadget3io/l50n288-phew-m5/"
# # fname = folder + "l50n288-phew-m5_108.starinfo"

# tab = genfromtxt(fname, names=True)
# primo = tab[tab['a_last'] == 0]

# msel = primo[primo['Mstar'] < 10.5]
# msel = msel[msel['Mstar'] > 10.0]
# # mratio = primo['WindMass'] / primo['Mass']
# mratio = msel['WindMass'] / msel['Mass']
# print sum(primo['WindMass']) / sum(primo['Mass'])
# print sum(msel['WindMass']) / sum(msel['Mass'])
# # da = primo['a_form'] - primo['a_acc']

# # step = 1000
# # plt.plot(da[::step], mratio[::step], "b.", alpha=0.2, markersize=6)
# # plt.axis([0.0, 0.1, 1.e-2, 3.0])
# # plt.yscale("log")
# # plt.show()

# ---------------- GASINFO ----------------
fig, ax = plt.subplots(1,1, figsize=(8,7))

# models = ['l25n144-phew-rcloud', 'l25n288-phew-m5']
# lgds = ['PhEW,25/144', 'PhEW,25/288']
models = ['l25n144-phew-rcloud', 'l25n144-phew-m5-spl', 'l25n288-phew-m5', 'l25n288-phew-m5-spl']
lgds = ['25/144', '25/144,Split', '25/288', "25/288,Split"]
ls = ["--", "-", "--", "-"]
clrs = ["blue", "blue", "red", "red"]

for mi, model in enumerate(models):
    folder = "/scratch/shuiyao/scidata/gadget3io/"+model+"/"
    fname = folder + model + "_098.gas.mh11"

    tab = genfromtxt(fname, names=True)
    primo = tab[tab['Mc'] == 0]

    if(SHOW_WIND_FRAC == True):    
        mratio = primo['WMass'] / primo['Mass']
    else:
        mratio = primo['Mass'] / Mgas_lres

    if(mi == 2 or mi == 3): mratio *= 8.0
    # da = primo['a_form'] - primo['a_acc']
    # print sum(primo['WMass']) / sum(primo['Mass'])

    hist, edges = histogram(mratio, bins=100)
    cen, cumfrac = reverse_cumsum(hist, edges)
    ax.plot(cen, cumfrac, clrs[mi], linestyle=ls[mi])

    # folder = "/scratch/shuiyao/scidata/gadget3io/"+model+"/"
    # fname = folder + model + "_108.starinfo.mh11"
    # tabsinfo = genfromtxt(fname, names=True)
    # primosfinfo = tabsinfo[tabsinfo['a_last'] == 0]
    # mratiosinfo = primosfinfo['WindMass'] / primosfinfo['Mass']
    # hist, edges = histogram(mratiosinfo, bins=100)
    # cen, cumfrac = reverse_cumsum(hist, edges)
    # ax.plot(cen, cumfrac, "r", linestyle=ls[mi])

    # latersfinfo = primosfinfo[primosfinfo['a_form'] > 0.5]
    # mratiosinfo = latersfinfo['WindMass'] / latersfinfo['Mass']
    # hist, edges = histogram(mratiosinfo, bins=100)
    # cen, cumfrac = reverse_cumsum(hist, edges)
    # ax.plot(cen, cumfrac, color="orange", linestyle=ls[mi])
    

lgd = legend.legend(ax)
lgd.loc = "upper right"
lgd.addLine((lgds[0], "blue", "--", 1))
lgd.addLine((lgds[1], "blue", "-", 1))
lgd.addLine((lgds[2], "red", "--", 1))
lgd.addLine((lgds[3], "red", "-", 1))
# lgd.addLine(("Halo Gas", "blue", "-", 1))
# lgd.addLine(("Stars", "red", "-", 1))
# lgd.addLine(("Stars (z<1)", "orange", "-", 1))
lgd.draw()

ax.set_ylim(0.0, 1.0)
if(SHOW_WIND_FRAC == True):
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel(r"$M_{w}/M_{p}$")
else:
    ax.set_xlim(0.0, 10.0)
    ax.set_xlabel(r"$M_{p}/M_{p,orig}$")
ax.set_ylabel("Fraction")
plt.title("z = 0.2")
plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
plt.show()

# plt.yscale("log")
# plt.show()

# ---------------- SFRINFO ----------------

# folder = "/proj/shuiyao/l25n288-phew-m5/SFRINFO/"
# fname = folder + "sfrinfo.0"

#acc = genfromtxt(fname, names=True)
# primo = acc[acc['LastSFTime'] == 0]
# merged = acc[acc['LastSFTime'] > 0]
# winds = acc[acc['LastSFTime'] < 0]

# mratio = primo['WindMass'] / primo['Mass']
# plt.hist(mratio, bins=linspace(0., 1., 100))
# plt.yscale("log")
# plt.show()

# print sum(primo['WindMass']) / sum(primo['Mass'])
# print sum(merged['WindMass']) / sum(merged['Mass'])
# print sum(winds['WindMass']) / sum(winds['Mass'])

print "Done."
