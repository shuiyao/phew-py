from mymod import *
from numpy import genfromtxt, isnan
import config_mpl

NBINS = 10
MMIN, MMAX = 11.0, 13.5
mcen = linspace(MMIN, MMAX, NBINS+1)
mcen = 0.5 * (mcen[:-1] + mcen[1:])

def find_mass_bin(mass, nbins=NBINS, mmin=MMIN, mmax=MMAX):
    dm = (mmax - mmin) / (float)(nbins)
    idx = (mass - mmin) / dm
    if(isnan(idx)): return -1
    if(idx < 0): idx = 0
    if(idx >= nbins-1): idx = nbins - 1
    return (int)(idx)

fbase = "/scratch/shuiyao/scidata/newwind/"
zstr = "z1"
models = ['l25n144-phew-m5kh30fs10', 'l25n144-phew-mach1', 'l25n144-phew-fa', 'l25n144-phew-norecouple']
lstyles = [':', '--', "--", "-"]
clrs = ['blue','orange','red',"magenta"]

fig, ax = plt.subplots(1,1, figsize=(6,6))
for mi in range(len(models)):
    if(not (mi == 2 or mi == 3)): continue
    flag_0 = array([0.0] * NBINS)
    flag_3 = array([0.0] * NBINS)
    flag_4 = array([0.0] * NBINS)
    Nnorm = array([0.0] * NBINS)
    fname = fbase + models[mi] + "/phewsinfo." + zstr 
    print "Reading file: ", fname
    tab = genfromtxt(fname, names=True)
    for part in tab:
        idx = find_mass_bin(part['LogMvir'])
        if(idx == -1): continue
        mass = part['Mass']
        Nnorm[idx] += mass
        flag = int(part['flag'])
        if(flag == 0 or flag == 1): flag_0[idx] += mass
        if(flag == 3): flag_3[idx] += mass
        if(flag == 4): flag_4[idx] += mass
    print "--- Model: ", models[mi], "---"
    print "Recouple: ", sum(flag_0) / sum(Nnorm) * 100.0, "%"
    print "Timeout: ", sum(flag_3) / sum(Nnorm) * 100.0, "%"
    print "Destroyed: ", sum(flag_4) / sum(Nnorm) * 100.0, "%"
    print "Stay: ", 100.0 - (sum(flag_0)+sum(flag_3)+sum(flag_4)) / sum(Nnorm) * 100.0, "%"
    flag_0 = flag_0 / Nnorm
    flag_3 = flag_3 / Nnorm
    flag_4 = flag_4 / Nnorm
    ax.plot(mcen, flag_0, color=clrs[0], linestyle=lstyles[mi])
    ax.plot(mcen, flag_3, color=clrs[1], linestyle=lstyles[mi])
    ax.plot(mcen, flag_4, color=clrs[2], linestyle=lstyles[mi])
    ax.plot(mcen, 1.0-(flag_0+flag_3+flag_4), color=clrs[3], linestyle=lstyles[mi])
ax.set_xlabel(r"$\log(M_\mathrm{vir}/M_\odot)$")
ax.set_ylabel("Fraction (Mass)")
ax.set_xlim(MMIN, MMAX)
ax.set_ylim(0.0, 1.0)
from pltastro import legend
lgd = legend.legend(ax)
lgd.fontsize = 12
lgd.addLine(("Recouple", clrs[0], "-", 1))
lgd.addLine(("Time out", clrs[1], "-", 1))
lgd.addLine(("Destroyed", clrs[2], "-", 1))
lgd.addLine(("Stay", clrs[3], "-", 1))
lgd.addLine((r"$\mathcal{M}_\mathrm{re}=1.0, f_\mathrm{re}=0.1$", "black", "--", 1))
lgd.addLine((r"$NoRec, f_\mathrm{re}=0.1$", "black", "-", 1))
lgd.loc = "upper left"
lgd.draw()
# lgd2 = legend.legend(ax)
# lgd2.fontsize = 12
# # lgd2.addLine((r"$\mathcal{M}_\mathrm{re}=2.0, f_\mathrm{re}=0.05$", "black", ":", 1))
# lgd2.addLine((r"$\mathcal{M}_\mathrm{re}=1.0, f_\mathrm{re}=0.1$", "black", "--", 1))
# lgd2.addLine((r"$NoRec, f_\mathrm{re}=0.1$", "black", "-", 1))
# lgd2.loc = "lower left"
# lgd2.draw()
plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
plt.show()
    
