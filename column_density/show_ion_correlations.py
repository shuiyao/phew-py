import ioformat
from mymod import *
import matplotlib.pyplot as plt

# x is perpendicular to the wind

fbase = "/scratch/shuiyao/Jneil/column_density/"
fnamex = fbase + "LowCond_v1700_chi300_cond_x.csv" # perpendicular
fnamey = fbase + "LowCond_v1700_chi300_cond.csv"

clrs_t = ["magenta", "red", "orange", "yellow"]

snapt = "t50"

if(snapt == "t90"): offset = 0
if(snapt == "t50"): offset = 20
if(snapt == "t25"): offset = 30

def convert_list(lst):
    lst = array(lst)
    lst = lst[lst > 0]
    lst = log10(lst)
    return lst

# 1   2   3  4   5      6    7    8     9    10
# OVI CIV NV CII NeVIII CIII MgII SiIII SiIV HI
ionnames = ["XX","OVI","CIV","NV","CII","NeVIII","CIII","MgII","SiIII","SiIV","HI"]
ionorder = [-1, 5, 3, -1, -1, 6, 2, 7, -1, 8, 0]

# cols = array([1, 5, 10]) # High-ions
# cols = array([2, 6, 7, 9, 10])

cols = array([1, 2, 5, 6, 7, 9, 10])
tab = ioformat.rcol(fnamey, cols + offset, separator=",", linestart=1)
for i in range(len(tab)): tab[i] = array(tab[i])
idxs = (tab[-1] > 1.e16)
for i in range(len(tab)): tab[i] = log10(tab[i][idxs])

# Strong Correlations
# HI: MgII, CIII       2
# OVI: NeVIII          1
# CIV: SiIV            1
# No/Weak Correlations
# HI: CIV, OVI         2
# CIV: CIII, OVI       2

from pltastro import frame, draw
frm = frame.multi(4,2)
pars = frm.params
pars.figsize = (8, 10)
pars.hspace = 0.40
pars.wspace = 0.25
pars.bottom = 0.08
pars.top = 0.95
pars.left = 0.1
pars.right = 0.9
pnls = frm.panels
pnls.xticksON = [True] * 8
pnls.yticksON = [True] * 8
fig, axs = draw(frm)

xions = [6, 6, 1, 0, 6, 6, 1, 1]
yions = [4, 3, 5, 2, 1, 0, 3, 0]
for i in range(8):
    xion, yion = xions[i], yions[i]
    axs[i].plot(tab[xion][::10], tab[yion][::10], "b.", markersize=4)
    axs[i].set_xlabel(ionnames[cols[xion]])
    axs[i].set_ylabel(ionnames[cols[yion]])    

plt.savefig("/scratch/shuiyao/figures/tmp.pdf")
plt.show()

# plt.show()

print "done"
