from mymod import *
import matplotlib.pyplot as plt

ionnames = ["XX","OVI","CIV","NV","CII","NeVIII","CIII","MgII","SiIII","SiIV","HI"]
fbase = "/scratch/shuiyao/Jneil/column_density/"
model = "LowCond_v1700_chi300_cond"
modelstr = "x300v1700"
# model = "HC_v1000_chi300_cond"
# modelstr = "x300v1000"
NPIX = 640000
fnamex = fbase + model + "_x.csv"
fnamey = fbase + model + ".csv"
ion = 10
cols = [20+ion]
ioncds = ioformat.rcol(fnamey, cols, separator=",", linestart=1)
ioncds = array(ioncds).reshape(800, 800)
ioncds = log10(ioncds)
plt.pcolor(ioncds, cmap=plt.get_cmap("jet"))
# plt.savefig("/scratch/shuiyao/figures/tmp.png")

print "done."
