# Find out info from the initwinds.*
from numpy import genfromtxt, linspace
import matplotlib.pyplot as plt
import config_mpl

fbase = "/proj/shuiyao/"
model = "l25n288-phew-m5"
# model = "l25n144-phew-rcloud"
unit_m = 3.4e6 # L25
NCPU = 256

models = ["l25n144-phew-rcloud", "l25n288-phew-m5"]
clrs = ["blue", "red"]

folder = fbase + model + "/WINDS/"
foutname = "windinit_"+model+".stat"

NBINS = 100.0

def find_abin_index(a, nbins=NBINS):
    if(a == 1.0): return 99
    return (int)(a / (1./NBINS))

unit_mass = 1.e10 / 0.7
def write_wind_stat():
    abins = linspace(1./NBINS, 1., NBINS)
    zbins = 1./abins - 1.
    mbins = [0.0] * (int)(NBINS)
    for icpu in range(NCPU):
        fname = folder + "initwinds." + str(icpu)
        print fname
        if(icpu == 0):
            tab = genfromtxt(fname, names=True)
            dtypes = tab.dtype
        else:
            tab = genfromtxt(fname, dtype=dtypes)
        for t in tab:
            mbins[find_abin_index(t['atime'])] += t['Mass'] * unit_mass

    fout = open(foutname, "w")
    fout.write("#atime z mass\n")
    for i in range((int)(NBINS)):
        fout.write("%7.5f %7.5f %7.5e\n" % (abins[i], zbins[i], mbins[i]))
    fout.close()

def show_wind_stat(fname, clr="black"):
    tab = genfromtxt(fname, names=True)
    p, = plt.plot(tab['atime'], tab['mass'], "-", color=clr)
    return p

lines = []    
for mi in range(len(models)):
    foutname = "windinit_"+models[mi]+".stat"
    p = show_wind_stat(foutname, clr=clrs[mi])
    lines.append(p)
plt.ylabel("Wind Mass")
plt.xlabel("ascale")
plt.legend(lines, ["25/144", "25/288"])
plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
plt.show()
    
print "done."
    
